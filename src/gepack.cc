/*
 *  Main authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Stefano Gualandi, 2012
 *
 *  Last update: October, 5h, 2012
 */

/// Boost Timer
#include <boost/progress.hpp>
using boost::timer;
timer  TIMER;

/// Limits on numbers
#include <limits>
using std::numeric_limits;

#include <vector>
using std::vector;

#include <stdio.h>
#include <fstream>

#include <gecode/driver.hh>
using namespace Gecode;
using namespace Gecode::Int;

#include "dag_pack.hh"

/// Main Script
class BinPacking : public Script {
   protected:
      /// 
      IntVarArray   x;  
      /// Load variable
      IntVarArray   l;
      /// Cost Variables
      IntVar        z;
   public:
      /// Actual model
      BinPacking( int n, int m, int C, const vector<int>& w, const vector< vector<int> >& D, int UB )
         :  x ( *this, n, 0, m-1 ), l ( *this, m, 0, C   ), z ( *this, 0, UB)
      {   
         /// Post binpacking constraint
         binpacking ( *this, l, x, w );
         /// Add costs to item-bin assignment
         IntVarArgs y(*this, n, 0, 1000000);
         for ( int i = 0; i < n; ++i ) {
            IntSharedArray S(m);
            for ( int j = 0; j < m; ++j )
               S[j] = D[i][j];
            element(*this, S, x[i], y[i]);
         }
         linear(*this, y, IRT_EQ, z);
         /// Make some random assignment
         //for ( int i = 0; i < n; ++i )
           // if ( w[i] > C/2 )
             //  rel( *this, x[i], IRT_EQ, x[i].min() );
      }

      BinPacking( bool share, BinPacking& s) : Script(share,s) {
         x.update ( *this, share, s.x );
      }

      /// Perform copying during cloning
      virtual Space*
         copy(bool share) {
            return new BinPacking(share, *this);
         }

      int totalDomain(void) {
         status();
         int dom = 0;
         for ( int i = x.size(); i > 0; --i )
            dom += x[i-1].size();
         fprintf(stdout,"DomTot %d (%d) %s ",dom, (x.size()*l.size()), (status() ? "OK" : "FAILED"));

         return status();
      }

      void getUBs( resources& U ) {
         for ( int i = 0; i < l.size(); ++i )
            U[i] = l[i].max();
      }

      void addArcs( DAG& G, int S, int T, int n, int m, int C, const vector<int>& w, 
            const vector< vector<int> >& D ) {
         for ( int i = 0; i < n-1; ++i ) 
            for ( IntVarValues j(x[i]); j(); ++j ) {
               for ( IntVarValues l(x[i+1]); l(); ++l ) {
                  if ( !(j.val() == l.val() && w[i] + w[i+1] > C) ) {
                     cost_t c = D[i+1][l.val()];
                     resources R(m,0);
                     R[l.val()] = w[i+1];
                     G.addArc( i*m+j.val(), (i+1)*m+l.val(), c, R );
                  }
               }
            }
         /// Arcs from the source node
         if ( w[0] <= C ) {
            for ( IntVarValues l(x[0]); l(); ++l ) {
               cost_t c = D[0][l.val()];
               resources R(m,0);
               R[l.val()] = w[0];
               G.addArc( S, l.val(), c, R );
            }
         }
         if ( w[n-1] <= C ) {
            for ( IntVarValues l(x[n-1]); l(); ++l ) {
               cost_t c = 0.0;
               resources R(m,0);
               G.addArc( (n-1)*m+l.val(), T, c, R );
               //fprintf(stdout,"(%d,%d)\t", (n-1)*m+l.val(), T);
            }
         }
      }
};



/// Find upper bounds to coloring
void singlePropagation ( int n, int m, int C, const vector<int>& w, const vector< vector<int> > D )
{
   timer TIMER;
   cost_t LB = 0;
   cost_t UB = 250;

   BinPacking* s = new BinPacking ( n, m, C, w, D, UB );
   
   if ( s->totalDomain() ) {
      /// My filtering 
      resources U(m, 0);
      s->getUBs(U); /// solo upper bounds
      int N = n*m+2;
      int M = n*m*m;
      int S = n*m;
      int T = n*m+1;
      DAG G (N, M, U);
      s->addArcs(G,S,T,n,m,C,w,D);
      /// Try to do filtering 
      if ( !G.isAcyclic() ) {
         fprintf(stdout, "It is not a DAG!\n");
         exit(EXIT_FAILURE);
      }


      fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
            LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

      G.filter(S,T,LB,UB);
      fprintf(stdout, "UB0 %.3f - ", UB);

      fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
            LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

      G.printArcs(n, m);
   }

}



/// MAIN PROGRAM

int main(int argc, char **argv)
{
   if ( argc != 2 )
      return EXIT_FAILURE;

   /// Beautify output
   std::ifstream infile(argv[1]); 
   if (!infile) 
   { 
      fprintf(stderr,"No %s file", argv[1]); 
      exit ( EXIT_FAILURE ); 
   }

   int n = 0;
   int m = 0;
   int C = 0;

   infile >> n >> C;
   vector<int> w(n,0);
   double w_tot = 0;
   for ( int i = 0; i < n; ++i ) {
      infile >> w[i];
      w_tot += w[i];
   }

   m = int(ceil(w_tot/C) + 5);

   vector< vector<int> > D;
   for ( int i = 0; i < n; ++i ) {
      vector<int> row;
      for ( int j = 0; j < m; j++ ) {
         row.push_back( rand() % 100 );
         fprintf(stdout,"%d ", row.back());
      }
      fprintf(stdout,"\n");
      D.push_back( row );
   }
      
   singlePropagation(n,m,C,w,D);

   fprintf(stdout,"%.3f\n", TIMER.elapsed());
   
   return EXIT_SUCCESS;
}

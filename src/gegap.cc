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

#include "cost_binpacking.hh"

#include <gecode/driver.hh>
using namespace Gecode;
using namespace Gecode::Int;

#include "dag_pack.hh"

/// Main Script
class BinPacking : public MinimizeScript {
   protected:
      /// 
      IntVarArray   x;  
      /// Load variable
      IntVarArray   l;
      /// Load variable
      IntVarArray   y;
      /// Cost Variables
      IntVar        z;
   public:
      /// Actual model
      BinPacking( int n, int m, int B, const vector<int>& w, const vector< vector<int> >& D )
         :  x ( *this, n, 0, m-1 ), 
            l ( *this, m, 0, B   ), 
            y ( *this, n, 0, Limits::max ), 
            z ( *this,    0, Limits::max )
      {   
         /// Post binpacking constraint
         binpacking ( *this, l, x, w );
         /// Add costs to item-bin assignment
         for ( int i = 0; i < n; ++i ) {
            IntSharedArray S(m);
            for ( int j = 0; j < m; ++j )
               S[j] = D[j][i];
            element(*this, S, x[i], y[i]);
         }
         linear(*this, y, IRT_EQ, z);
         /// Post cost bin packing constraints
         //IntSharedArray C(n*m);
         //int idx = 0;
         //for ( int i = 0; i < n; ++i ) {
            //for ( int j = 0; j < m; ++j, ++idx ) {
               //C[i*m+j] = D[j][i];
               //fprintf(stdout,"%d %d %d %d %d\t", j, i, idx, C[i*m+j], C[idx]);
            //}
            //fprintf(stdout,"\n");
         //}
         //cost_binpacking(*this, l, x, w, C);

         branch(*this, x, INT_VAR_SIZE_DEGREE_MAX, INT_VAL_MIN);
      }
      /// Constructor for cloning \a s
      BinPacking( bool share, BinPacking& s) : MinimizeScript(share,s) {
         x.update ( *this, share, s.x );
         l.update ( *this, share, s.l );
         y.update ( *this, share, s.y );
         z.update ( *this, share, s.z );
      }
      /// Perform copying during cloning
      virtual Space* copy(bool share) {
         return new BinPacking(share, *this);
      }
      /// Return cost
      virtual IntVar cost(void) const {
         return z;
      }
  
      virtual int getValueSL(void) const { return z.val(); }
      virtual int getValueUB(void) const { return z.max(); }
      virtual int getValueLB(void) const { return z.min(); }

      int totalDomain(void) {
         int dom = 0;
         for ( int i = 0; i < x.size(); ++i ) 
            for ( IntVarValues v(x[i]); v(); ++v, dom++ );
         fprintf(stdout,"DomTot %d %d#%d %s ", dom, z.min(), z.max(), (status() ? "OK" : "FAILED"));
         return status();
      }

      void getUBs( resources& U ) {
         for ( int i = 0; i < l.size(); ++i )
            U[i] = resource_t(l[i].max());
      }

      void getLBs( resources& L ) {
         for ( int i = 0; i < l.size(); ++i )
            L[i] = resource_t(l[i].min());
      }

      void addArcs( DAG& G, int S, int T, int n, int m, int C, 
            const vector<int>& w, 
            const vector< vector<int> >& D ) {
         for ( int i = 0; i < n-1; ++i ) 
            for ( IntVarValues j(x[i]); j(); ++j ) {
               for ( IntVarValues l(x[i+1]); l(); ++l ) {
                  if ( !(j.val() == l.val() && w[i] + w[i+1] > C) ) {
                     cost_t c = D[l.val()][i+1];
                     resources R(m,0);
                     R[l.val()] = x[i+1].size();
                     G.addArc( i*m+j.val(), (i+1)*m+l.val(), c, R );
                  }
               }
            }
         /// Arcs from the source node
         for ( IntVarValues h(x[0]); h(); ++h ) 
            if ( w[0] <= l[h.val()].max() ) {
               cost_t c = D[h.val()][0];
               resources R(m,0);
               R[h.val()] = w[0];
               G.addArc( S, h.val(), c, R );
            }

         for ( IntVarValues h(x[n-1]); h(); ++h ) 
            if ( w[n-1] <= l[h.val()].max() ) {
               cost_t c = 0.0;
               resources R(m,0);
               G.addArc( (n-1)*m+h.val(), T, c, R );
            }
         
      }
};

/// SOlve the instance only using CP (no reduced-cost filtering)
void onlyCP ( int n, int m, int C, const vector<int>& w, const vector< vector<int> > D )
{
  /// Solution of the problem
  Support::Timer t;
  t.start();

  Search::Options so;
  
  BinPacking* s = new BinPacking ( n, m, C, w, D );

  BAB<BinPacking> e(s, so);
  int UB = -1;
  int tic = 0;
  do {
     BinPacking* ex = e.next();
     if ( ex == NULL )
        break;
     UB = ex->getValueSL();
     if ( tic++ % 13 == 0 )
        fprintf(stdout,"UB %d Nodes %ld Memory %ld Time %.3f\n", 
              UB, e.statistics().node, e.statistics().memory, TIMER.elapsed());
     delete ex;
  } while(true);

  fprintf(stdout,"SL %d Nodes %ld Memory %ld Time %.3f\n", 
        UB, e.statistics().node, e.statistics().memory, TIMER.elapsed());
}

/// Find upper bounds to coloring
void singlePropagation ( int n, int m, int C, const vector<int>& w, const vector< vector<int> > D )
{
   timer TIMER;

   BinPacking* s = new BinPacking ( n, m, C, w, D );
   
   if ( s->totalDomain() ) {
      cost_t LB = s->getValueLB();
      cost_t UB = s->getValueUB();
      /// My filtering 
      resources U(m, 0);
      resources L(m, 0);
      s->status();
      s->getUBs(U); /// solo upper bounds
      s->getLBs(L);
      node_t N = n*m+2;
      edge_t M = n*m*m;
      node_t S = n*m;
      node_t T = n*m+1;
      DAG G (N, M, U);
      s->addArcs(G,S,T,n,m,C,w,D);

      fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
            LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

      G.filter(S,T,LB,UB);
      G.printArcs(n, m);

      fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
            LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);
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
   //# FORMAT OF INSTANCES:
   //#
   //# resource capacity of agent i (i=1,...,m)
   
   int n = 0;  /// items
   int m = 0;  /// bins (?)
   int C = 0;

   // number of agents (m), number of jobs (n)
   infile >> m >> n;
   // for each agent i (i=1,...,m) in turn:
   //   cost of allocating job j to agent i (j=1,...,n)
   vector< vector<int> > D;
   for ( int j = 0; j < m; ++j ) {
      vector<int> row(n,0);
      for ( int i = 0; i < n; ++i ) 
         infile >> row[i];
      D.push_back( row );
   }
   // for each agent i (i=1,...,m) in turn:
   //   resource consumed in allocating job j to agent i (j=1,...,n)
   vector<int> w(n,0);
   double w_tot = 0.0;
   for ( int i = 0; i < n; ++i ) {
      infile >> w[i];
      w_tot += w[i];
   }
   int tmp;
   for ( int j = 0; j < m-1; ++j ) 
      for ( int i = 0; i < n; ++i )
         infile >> tmp;
   infile >> C;
   for ( int j = 0; j < m-1; ++j ) 
      infile >> tmp; 

   //for ( int j = 0; j < m; ++j ) {
      //for ( int i = 0; i < n; ++i )
         //fprintf(stdout,"%d ",D[j][i]);
      //fprintf(stdout, "\n");
   //}

   int C0 = C;
   while ( w_tot/m > C )
      C++;

   fprintf(stdout,"%d %d %d %d\n", n, m, C0, C);
   onlyCP(n,m,C,w,D);
   //singlePropagation(n,m,C,w,D);

   fprintf(stdout,"%.3f\n", TIMER.elapsed());
   
   return EXIT_SUCCESS;
}

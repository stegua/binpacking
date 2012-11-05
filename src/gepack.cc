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

/// Cutoff object for the script
class MyCutoff : public Search::Stop {
   private:
      Search::TimeStop*    ts; 
      Search::MemoryStop*  ms; 

      MyCutoff(unsigned int time, unsigned int memory)
         : ts((time > 0)   ? new Search::TimeStop(time) : NULL),
         ms((memory > 0) ? new Search::MemoryStop(memory) : NULL) {}
   public:
      virtual bool stop(const Search::Statistics& s, const Search::Options& o) {
         return
            ((ts != NULL) && ts->stop(s,o)) ||
            ((ms != NULL) && ms->stop(s,o));
      }
      virtual bool stopTime(const Search::Statistics& s, const Search::Options& o) {
         return
            ((ts != NULL) && ts->stop(s,o));
      }
      virtual bool stopMemory(const Search::Statistics& s, const Search::Options& o) {
         return
            ((ms != NULL) && ms->stop(s,o));
      }
      static Search::Stop*
         create(unsigned int time, unsigned int memory) {
            if ((time == 0) && (memory == 0))
               return NULL;
            else
               return new MyCutoff(time,memory);
         }
      ~MyCutoff(void) {
         delete ts; delete ms;
      }
};

/// Main Script
class BinPacking : public Script {
   protected:
      /// 
      IntVarArray   x;  
      /// Load variable
      IntVarArray   l;
   public:
      /// Actual model
      BinPacking( int n, int m, int C, const vector<int>& w )
         :  x ( *this, n, 0, m-1 ), l ( *this, m, 0, C   )
      {   
         /// Post binpacking constraint
         binpacking ( *this, l, x, w );
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

      void addArcs( DAG& G, int S, int T, int n, int m, int C, const vector<int>& w ) {
         for ( int i = 0; i < n-2; ++i ) 
            if ( w[i] + w[i+1] <= C) {
               for ( IntVarValues j(x[i]); j(); ++j ) {
                  for ( IntVarValues l(x[i+1]); l(); ++l ) {
                     if ( j.val() != l.val() ) {
                        cost_t c = 1.0;
                        resources R(m,0);
                        R[l.val()] = w[i+1];
                        G.addArc( i*m+j.val(), (i+1)*m+l.val(), c, R );
                     }
                  }
               }
            }
         /// Arcs from the source node
         if ( w[0] <= C ) {
            for ( IntVarValues l(x[0]); l(); ++l ) {
               cost_t c = 1.0;
               resources R(m,0);
               R[l.val()] = w[0];
               G.addArc( S, l.val(), c, R );
            }
         }
         if ( w[n-1] <= C ) {
            for ( IntVarValues l(x[n-1]); l(); ++l ) {
               cost_t c = 1.0;
               resources R(m,0);
               G.addArc( (n-1)*m+l.val(), T, c, R );
               //fprintf(stdout,"(%d,%d)\t", (n-1)*m+l.val(), T);
            }
         }
      }
};



/// Find upper bounds to coloring
void singlePropagation ( int n, int m, int C, const vector<int>& w )
{
   BinPacking* s = new BinPacking ( n, m, C, w );
   if ( s->totalDomain() ) {
      /// My filtering 
      resources U(m, 0);
      s->getUBs(U); /// solo upper bounds
      int N = n*m+2;
      int M = n*m*m;
      int S = n*m;
      int T = n*m+1;
      DAG G (N, M, U);
      s->addArcs(G,S,T,n,m,C,w);
      /// Try to do filtering 
      if ( !G.isAcyclic() ) {
         fprintf(stdout, "It is not a DAG!\n");
         exit(EXIT_FAILURE);
      }

      timer TIMER;
      cost_t LB = 0;
      cost_t UB = n+2;

      fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
            LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

      //G.filter(S,T,LB,UB);
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
   int w_tot = 0;
   for ( int i = 0; i < n; ++i ) {
      infile >> w[i];
      w_tot += w[i];
   }

   m = w_tot/C + 2;

   singlePropagation(n,m,C,w);

   fprintf(stdout,"%.3f\n", TIMER.elapsed());
   
   return EXIT_SUCCESS;
}

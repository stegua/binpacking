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
      BinPacking( int n, int m, int B, const vector<int>& w, const vector< vector<int> >& D, bool cost )
         :  x ( *this, n, 0, m-1 ), 
            l ( *this, m, 0, B   ), 
            y ( *this, n, 0, Limits::max ), 
            z ( *this,    0, 2000)//Limits::max )
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
         IntSharedArray C(n*m);
         for ( int i = 0; i < n; ++i ) 
            for ( int j = 0; j < m; ++j ) 
               C[i*m+j] = D[j][i];
         /// If to use the new constraint 
         if ( cost )
            cost_binpacking(*this, l, x, z, w, C);

         //IntVarArgs zv(1);
         //zv[0] = z;
         //branch(*this, zv, INT_VAR_SIZE_MIN, INT_VAL_MED);
        
         branch(*this, x,  INT_VAR_SIZE_DEGREE_MAX, INT_VAL_MIN);
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
};

/// SOlve the instance only using CP (no reduced-cost filtering)
void onlyCP ( int n, int m, int C, const vector<int>& w, const vector< vector<int> > D, bool cost=true )
{
  /// Solution of the problem
  Support::Timer t;
  t.start();

  Search::Options so;
  
  BinPacking* s = new BinPacking ( n, m, C, w, D, cost );

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

/// MAIN PROGRAM

int main(int argc, char **argv)
{
   if ( argc < 2 )
      return EXIT_FAILURE;
   
   bool cost = !(argc == 3);

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
   onlyCP(n,m,C,w,D,cost);

   fprintf(stdout,"%.3f\n", TIMER.elapsed());
   
   return EXIT_SUCCESS;
}

/*
 *  Main authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Stefano Gualandi, 2013
 *
 *  Last update: April, 27th, 2013
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

#include "multibin-packing.hh"

#include <gecode/driver.hh>
using namespace Gecode;
using namespace Gecode::Int;


/// Main Script
class MultiBinPacking : public Script {
   protected:
      /// x_i = j : item i in bin j
      IntVarArray   x;  
      /// Load variable
      IntVarArray   y;
   public:
      /// Actual model
      MultiBinPacking( int n, int m, int k, const vector<int>& b, const vector< vector<int> >& A )
         :  x ( *this, n,   0, m-1         ), 
            y ( *this, m*k, 0, Limits::max ) 
      { 
         /// Every bin has the same capacity for a given dimension
         IntSharedArray B(k);
         for ( int l = 0; l < k; ++l )
            B[l] = b[l];

         /// Item size for each dimension
         IntSharedArray C(n*k);
         for ( int i = 0; i < n; ++i ) 
            for ( int l = 0; l < k; ++l ) 
               C[i*k+l] = A[i][l];

         /// Post multidimensional binpacking constraint
         multibinpacking(*this, n, m, k, y, x, C, B);

         /// Branching strategy 
         branch(*this, x, INT_VAR_DEGREE_SIZE_MAX(), INT_VAL_MIN());
      }
      /// Constructor for cloning \a s
      MultiBinPacking( bool share, MultiBinPacking& s) : Script(share,s) {
         x.update ( *this, share, s.x );
         y.update ( *this, share, s.y );
      }
      /// Perform copying during cloning
      virtual Space* copy(bool share) {
         return new MultiBinPacking(share, *this);
      }
      void printSol() {
         int n = x.size();
         for ( int i = 0; i < n; ++i )
            fprintf(stdout, "%d,", x[i].val());
         fprintf(stdout,"\n");
      }
};

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

/// Solve the instance only using CP (no reduced-cost filtering)
void onlyCP ( int n, int m, int k, const vector<int>& b, const vector< vector<int> >& A )
{
  /// Solution of the problem
  Support::Timer t;
  t.start();

  Search::Options so;
 
  so.c_d = 1;
  so.threads = 1;

  MultiBinPacking* s = new MultiBinPacking ( n, m, k, b, A );

  so.stop = MyCutoff::create( 120*1000, 1e+09 );
  
  DFS<MultiBinPacking> e(s, so);
     
  MultiBinPacking* ex = e.next();
  if ( e.stopped() ) {
     fprintf(stdout, "%.1f \tWARNING: STOPPED, IT IS ONLY AN UPPER BOUND! ",TIMER.elapsed());
     MyCutoff* myc = dynamic_cast<MyCutoff *>(so.stop);
     if ( myc->stopTime(e.statistics(),so) )
        fprintf(stdout," TIME LIMIT!\n");
     if ( myc->stopMemory(e.statistics(),so) )
        fprintf(stdout," MEMORY LIMIT!\n");
  }
 
  if ( ex == NULL )
     printf("NO SOLUTION\t");
  else
     ex->printSol();
  delete ex;

  fprintf(stdout,"Nodes %ld Memory %ld Time %.3f\n", 
        e.statistics().node, e.statistics().memory, TIMER.elapsed());
}

/// MAIN PROGRAM

int main(int argc, char **argv)
{
   if ( argc < 2 )
      return EXIT_FAILURE;

   /// Beautify output
   std::ifstream infile(argv[1]); 
   if (!infile) 
   { 
      fprintf(stderr,"No %s file", argv[1]); 
      exit ( EXIT_FAILURE ); 
   }
   
   int n = 0;  /// items
   int m = 0;  /// bins 
   int k = 0;  /// dimension
   infile >> m >> n >> k;

   vector< vector<int> > A;
   for ( int i = 0; i < n; ++i ) {
      vector<int> row(k,0);
      for ( int l = 0; l < k; ++l ) 
         infile >> row[l];
      A.push_back( row );
   }
   
   vector<int> b(k,0);
   for ( int l = 0; l < k; ++l ) 
      infile >> b[l];
   
   onlyCP(n,m,k,b,A);

   return EXIT_SUCCESS;
}

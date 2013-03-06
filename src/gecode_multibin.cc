/*
 *  Main authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Stefano Gualandi, 2013
 *
 *  Last update: March, 7t, 2013
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
#include <gecode/minimodel.hh>
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
         Matrix<IntVarArray> L(y, m, k);
         /// Capacity constraint for each dimension
         for ( int l = 0; l < k; ++l )
            rel ( *this, L.row(l), IRT_LQ, b[l] );
         /// Post binpacking constraints
         for ( int l = 0; l < k; ++l ) {
            IntArgs s(n);
            for ( int i = 0; i < n; ++i )
               s[i] = A[i][l];
            binpacking ( *this, L.row(l), x, s );
         }
         if ( status() == SS_FAILED )
            printf("FAILED ROOT\n");
         for ( int j = 0; j < m; j++ ) {
            for ( int l = 0; l < k; ++l )
               printf("%d#[%d,%d]  ", b[l], L(j,l).min(), L(j,l).max());
            printf("\n");
         }
         /// Branching strategy 
//         branch(*this, x, INT_VAR_SIZE_DEGREE_MAX, INT_VAL_MIN);
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
  
  MultiBinPacking* s = new MultiBinPacking ( n, m, k, b, A );

  so.stop = MyCutoff::create( 3600*1000, 1e+09 );
  
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
   
   vector<int> b(m,0);
   for ( int j = 0; j < m; ++j ) 
      infile >> b[j];
   
   onlyCP(n,m,k,b,A);

   fprintf(stdout,"%.3f\n", TIMER.elapsed());
   
   return EXIT_SUCCESS;
}

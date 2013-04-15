/// My typedefs
#include "dag_pack.hh"

/// Boost Timer
#include <boost/progress.hpp>
using boost::timer;

using namespace boost;

#include <fstream>

using namespace boost;

///------------------------------------------------------------------------------------------
/// Main function
int
main (int argc, char **argv)
{
   if ( argc != 2 ) {
      fprintf(stdout, "usage: ./multibin <filename>\n");
      exit ( EXIT_FAILURE );
   }

   /// Read instance from the OR-lib
   std::ifstream infile(argv[1]); 
   if (!infile) 
      exit ( EXIT_FAILURE ); 
  
   /// First read the data instace
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
  
   /// ... then build the corresponding graph
   int N = n*m+2;
   int M = (n-1)*(m*m)+2*m;
   int K = k*m;
   int S = N-2;
   int T = N-1;

   resources U(K,0);
   for ( int j = 0; j < m; ++j )
      for ( int l = 0; l < k; ++l )
         U[j*k+l] = b[l];

   fprintf(stdout,"\n");
   for ( int l = 0; l < K; ++l )
      fprintf(stdout,"%d ", U[l]);
   fprintf(stdout,"\n");

   DAG G (N, M, U);
   /// Archi da item a item
   for ( int i = 0; i < n-1; ++i ) {
      for ( int j = 0; j < m; ++j ) {
         resources R(K,0);
         for ( int l = 0; l < k; ++l )
            R[j*k+l] = A[i+1][l];
         for ( int h = 0; h < m; ++h ) {
            //for ( int ll = 0; ll < K; ++ll )
               //fprintf(stdout,"%d ", R[ll]);
            //fprintf(stdout,"\t (%d, %d)\n", i*m+h, ((i+1)*m)+j);
            if ( j == 0 )
               G.addArc( i*m+h, ((i+1)*m)+j, A[i+1][0], R);
            else
               G.addArc( i*m+h, ((i+1)*m)+j, 0, R);
         }
      }
   }      
   
   /// Archi dalla sorgente 
   for ( int j = 0; j < m; ++j ) {
      resources R(K,0);
      for ( int l = 0; l < k; ++l )
         R[j*k+l] = A[0][l];
      //for ( int ll = 0; ll < K; ++ll )
         //fprintf(stdout,"%d ", R[ll]);
      //fprintf(stdout,"\t (%d, %d)\n", S, j);
      if ( j == 0 )
         G.addArc( S, j, A[0][0], R);
      else
         G.addArc( S, j, 0, R);
   }
   /// Archi dalla destinazione
   for ( int j = 0; j < m; ++j ) {
      resources R(K,0);
      //fprintf(stdout,"(%d, %d)\n", (n-1)*m+j, T);
      G.addArc( (n-1)*m+j, T, 0, R);
   }

   if ( !G.isAcyclic() ) {
      fprintf(stdout, "It is not a DAG!\n");
      exit(EXIT_FAILURE);
   }

   timer TIMER;
   cost_t LB = 0;
   cost_t UB = b[0];
   
   //fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
     //    LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

   int status = G.filter(S,T,LB,UB);
   //fprintf(stdout, "UB0 %.3f - ", UB);
   
   fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d) - Status %d\n", 
         LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N, status);

   //G.printArcs(n, m);

   return 0;
}

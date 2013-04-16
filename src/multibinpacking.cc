/// My typedefs
#include "dag_pack.hh"

/// Boost Timer
#include <boost/progress.hpp>
using boost::timer;

using namespace boost;

#include <fstream>

using namespace boost;

#include <map>
using std::map;

///------------------------------------------------------------------------------------------
/// Main function
int
main (int argc, char **argv)
{
   if ( argc == 1 ) {
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

   //fprintf(stdout,"\n");
   //for ( int l = 0; l < K; ++l )
      //fprintf(stdout,"%d ", U[l]);
   //fprintf(stdout,"\n");

   DAG G (N, M, U);

   typedef std::pair<int, int> NodePair;
   map< NodePair, edge_t > As;
   /// Archi da item a item
   for ( int i = 0; i < n-1; ++i ) {
      for ( int j = 0; j < m; ++j ) {
         resources R(K,0);
         for ( int l = 0; l < k; ++l )
            R[j*k+l] = A[i+1][l];
         for ( int h = 0; h < m; ++h ) {
            //if ( j == 0 )
               //G.addArc( i*m+h, ((i+1)*m)+j, A[i+1][0], R);
            //else
            edge_t arc_id = G.addArc( i*m+h, ((i+1)*m)+j, 0, R);
            As[make_pair(i*m+h, ((i+1)*m+j))] = arc_id;
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
      //if ( j == 0 )
         //G.addArc( S, j, A[0][0], R);
      //else
         edge_t arc_id = G.addArc( S, j, 0, R);
         As[make_pair(S, j)] = arc_id;
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
   
   cost_t LB;
   cost_t UB;
   int status;
   int ris = 0;
   /// Cicla sulle risorse
   //int l = atoi(argv[2]);
   for ( int q = 0; q < 1; ++q ) {
      for ( int l = 0; l < k; ++l ) {
         ris++;
         fprintf(stdout, "\n");
         LB = 0;
         UB = b[l];
   
         /// Archi da item a item
         for ( int i = 0; i < n-1; ++i ) {
            for ( int j = 0; j < m; ++j ) {
               for ( int h = 0; h < m; ++h ) {
                  if ( j == q )
                     G.setArcCost(As[make_pair(i*m+h, ((i+1)*m+j))], A[i+1][l]);
                  else
                     G.setArcCost(As[make_pair(i*m+h, ((i+1)*m+j))], 0);
               }
            }
         }      

         /// Archi dalla sorgente 
         for ( int j = 0; j < m; ++j ) {
            if ( j == q )
               G.setArcCost( As[make_pair(S, j)], A[0][l]);
            else
               G.setArcCost( As[make_pair(S, j)], 0);
         }

         status = G.filter(S,T,LB,UB);

         fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d) - Status %d\n", 
               LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N, status);
         if ( status == 1 )
            goto QUIT;
      }
   }
QUIT:
   fprintf(stdout,"LOG \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d) - Status %d - Ris %d\n", 
         TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N, status, ris);

   return 0;
}

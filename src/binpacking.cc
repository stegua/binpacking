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
      fprintf(stdout, "usage: ./rcsp <filename>\n");
      exit ( EXIT_FAILURE );
   }

   /// Read instance from the OR-lib
   std::ifstream infile(argv[1]); 
   if (!infile) 
      exit ( EXIT_FAILURE ); 
   
   //int n = 5; /// Number of items
   //int m = 2; /// Number of bins
   //int C = 7; /// Bin capacity
   //int w[] = {4,4,2,2,2}; /// Items' weights
   int n = 10; /// Number of items
   int m = 3; /// Number of bins
   //int C = 100; /// Bin capacity
   //int D = 99;
   //int w[] = {49,41,34,33,29,26,26,22,20,19};

   int M, S, T;
   infile >> n >> M >> m >> S >> T;
   int k = 2*m;
   /// Read lower and upper bund resource consumptions
   resources U(k,0);
   for ( int l = 0; l < k; ++l )
      infile >> U[l];
   
   /// Build the graph (the beasley instances have only upper bounds)
   int N = n*m+2;
   DAG G (N, M, U);

   int v;
   int w;
   int c;
   for ( int e = 0; e < M; ++e ) {
      infile >> v >> w >> c;
      resources R(k,0);
      for ( int l = 0; l < k; ++l )
         infile >> R[l];
      G.addArc(v, w, 1, R);
   }

   //for ( int i = 0; i < n-1; i++ ) {
      //for ( int j = 0; j < m; j++ ) {
         //for ( int l = 0; l < k; ++l ) {
            //resources Ri(2*k,0);
            //Ri[l] = w[i+1];
            //Ri[l+m] = -w[i+1];
            //G.addArc((i*m)+j, ((i+1)*m)+l, 1, Ri);
            //fprintf(stdout,"%d %d %.f\n", i*m+j, (i+1)*m+l,Ri[l]);
         //}
      //}
      //fprintf(stdout,"\n");
   //}
//
   //for ( int j = 0; j < m; ++j ) {
    //if ( j == 0 ) {
       //resources Rs(k,0);
      //Rs[j] = w[0];
      //Rs[j+m] = -w[0];
      //G.addArc(S,j,1,Rs);
      //fprintf(stdout,"%d %d\t", S, j);
    //} 
      //resources Rt(k,0);
      //G.addArc((n-1)*m+j, T, 1, Rt);
      //fprintf(stdout,"%d %d\t", (n-1)*m+j, T);
      //
      //fprintf(stdout,"\n");
   //}

   if ( !G.isAcyclic() ) {
      fprintf(stdout, "It is not a DAG!\n");
      exit(EXIT_FAILURE);
   }

   timer TIMER;
   cost_t LB = 0;
   cost_t UB = n+2;
   
   fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
         LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

   G.filter(S,T,LB,UB);
   fprintf(stdout, "UB0 %.3f - ", UB);
   
   fprintf(stdout,"LB %.3f \t UB %.3f \t Time %.3f \t Arc removed %d (%d) - Nodes %d (%d)\n", 
         LB, UB, TIMER.elapsed(), int(G.arcsLeft()), int(G.num_arcs()), G.num_nodes(), N);

   G.printArcs(n, m);

   return 0;
}

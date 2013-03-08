#include "CBack.h"

#include "graph.h"
#include "reorder.h"

/// Complementary graph 
int* Gb;
int* Eb;

int n, LB;
int* S;
int* T; /// table of max clique sizes (one for each vertex)

void PrintCount()
{ 
   int i;
//   for ( i = 0; i < n; i++ )
 //     printf ("T[%d]=%d\n",i,T[i]);
   printf("\nw(G)=%d.\tClique: \t",LB); 
   for ( i = 0; i < n; i++ )
      if ( S[i] == 1 )
         printf("%d ", i);
   printf("\n");
   free(S);
   free(T);
   free(Eb);
   free(Gb);
}

void set_graph(int argc, char* argv[]) {
   int p;
   int i, j, m, mb;
	FILE *fp;
	graph_t *g;
   int* table;

   if ( argc != 2 )
      exit(2);

   fp=fopen(argv[1],"rb");
   if (fp==NULL) {
      perror(argv[1]);
      exit(2);
   }
   g=graph_read_dimacs(fp);
	if (g==NULL) {
		fprintf(stderr,"Error in graph file.\n");
      exit(2);
	}
	fclose(fp);

   table=reorder_by_greedy_coloring(g,FALSE);

   /// Metti il grafo nel formato che voglio
   n = g->n;   
   m = graph_edge_count(g);
   mb = ((n-1)*n/2-m);
   S = (int*) malloc(n*sizeof(int));
   T = (int*) malloc(n*sizeof(int));
   
   Eb = (int*) malloc((n+1)*sizeof(int));
   Gb = (int*) malloc(( mb)*sizeof(int));

   Eb[0]=-1;
	for (i=0; i < n; i++) {
      S[i]=0;
      T[i]=0;
		for (j=i+1; j < n; j++) {
//			if ( !SET_CONTAINS_FAST(g->edges[table[i]],table[j]))
			if ( !SET_CONTAINS_FAST(g->edges[table[n-1-i]],table[n-1-j]))
            Gb[p++]=j;
      }
      Eb[i+1]=p-1;
   }

   LB=1;
   
   graph_free(g);
   free(table);
}

void MaxClique(int argc, char* argv[]) {
   int* C;
   int v, w, c, s, out_clique, l;
   
   set_graph(argc, argv);
   Fiasco = PrintCount;
   
   C = (int*) Ncalloc(n, sizeof(int));
   for ( v = 0; v < n; v++ )
      C[v] = 0;
   l = n-Choice(n);
   C[l] = 1;
   T[l] = LB;
   s=1;
   out_clique=l;
   for ( w = Eb[l+1]; w > Eb[l] ; w-- ) {
      C[Gb[w]] = 2;
      out_clique--;
   }
   for (v = l+1; v < n; v++ ) {
      /// Skip removed vertices
      if ( C[v] < 2 && s+T[v] > LB ) {
         if ( s+1+T[v] <= LB )
            Backtrack();
         if ( s + n-out_clique <= LB  )
            Backtrack();
         c = 2 - Choice(2);
         C[v] = c;
         if ( c == 1 ) {
            /// Remove neighbors from candidate list
            for ( w = Eb[v+1]; w > Eb[v] ; w-- ) 
               if ( C[Gb[w]] == 0 ) {
                  out_clique++;
                  C[Gb[w]] = 2;
               }
            /// Update clique size and clique vector
            s++;
            if ( s > LB ) {
               printf ("Node: %d | w(G)=%d | \n",l,s);
               LB = s;
               for ( w = l; w < n; w++ )
                  S[w] = C[w];
            }
         } else
            out_clique++;
      }
   }
 //  printf ("Node: %d | w(G)=%d | \n",l,s);
   T[l] = LB;
   Backtrack();
}

int main(int argc, char* argv[]) Backtracking(MaxClique(argc,argv))

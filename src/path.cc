#include "path.hh"
#include "dag_pack.hh"

Path::Path ( const DAG& G, node_t Source, node_t Target, const vector< edge_t >& P ) {
   int k = G.getKap();
   Rc.resize(k,0);
   /// Set the path cost
   c = 0; 
   d = 0;
   /// Insert the nodes in the path
   node_t pred = (G.getArc(P[Target])).w;
   while ( pred != Source ) {
      path.push_back(P[pred]);
      c += (G.getArc(P[pred])).c;
      d += (G.getArc(P[pred])).d;
      for ( int l = 0; l < k; ++l )
         Rc[l] += (G.getArc(P[pred])).r[l];
      pred = (G.getArc(P[pred])).v;
   }
}

cost_t
Path::updateCost( const DAG& G ) {
   int np = path.size();
   c = 0;
   for ( int i = 0; i < np; ++i ) 
      c += (G.getArc(path[i])).c;
   return c; 
}


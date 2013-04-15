#include "path.hh"

Path::Path ( const DAG& G, node_t Source, node_t Target, const vector< edge_t >& P ) {
   int k = G.getKap();
   Rc.resize(k,0);
   /// Set the path cost
   c = 0; 
   d = 0;
   /// Insert the nodes in the path
   node_t pred = (G.getArc(P[Target])).w;
   while ( pred != Source ) {
      path.push_front(pred);
      c += (G.getArc(P[pred])).c;
      d += (G.getArc(P[pred])).d;
      for ( int l = 0; l < k; ++l )
         Rc[l] += (G.getArc(P[pred])).r[l];
      pred = (G.getArc(P[pred])).v;
   }
   path.push_front(Source);
}

/*Path::Path ( node_t Source, node_t Target, cost_t _c, const Arc& a, const vector< node_t >& Pf, const vector< node_t >& Pb ) {
   /// Set the path cost
   c = _c; 
   /// Insert the nodes on the arcs
   path.push_front(a.w);
   path.push_front(a.v);
   /// Given the arc a=(v,w) from v go back to the source ...
   {
      node_t cur = a.v;
      while ( cur != Source ) {
         node_t pre = Pf[cur];
         path.push_front(pre);
         cur = pre;
      }
   }
   /// ... and from go forward to the target
   {
      node_t cur = a.w;
      while ( cur != Target ) {
         node_t suc = Pb[cur];
         path.push_back(suc);
         cur = suc;
      }
   }
}

Path::Path ( node_t Target, cost_t _c, const stack<node_t>& theStack ) {
   /// Set the path cost
   c = _c; 
   /// Insert the nodes in the path
   stack<node_t> MyStack(theStack);
   path.push_front(Target);
   while(!MyStack.empty()) {
      path.push_front(MyStack.top());
      MyStack.pop();
   }
}*/


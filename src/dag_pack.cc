#include "dag_pack.hh"
#include "cost_resources_pack.hh"
#include <cmath>

/// Remove all the arcs incidient in the node given by <vit>
void
DAG::clearVertex(NodeIter vit) {
   for ( FSArcIterPair it = vit->getIterFS(); it.first != it.second; ) {
      FSArcIter a = it.first;
      ++it.first;  /// First increment the pointer, for safe in-place removals
      removeOppositeArc(a);
   }
   for ( BSArcIterPair it = vit->getIterBS(); it.first != it.second; ) {
      BSArcIter a = it.first;
      ++it.first;  /// First increment the pointer, for safe in-place removals
      removeOppositeArc(a);
   }
   N.erase(vit);
}

/// Run topological sort and return a vector of topological sorted nodes
/// [Note: we can do also topological sort via DFS; this can be useful to
///  explore the verteces in random order... maybe we get a different path
///  at eah iteration: give it a try!]
vector<node_t>
DAG::topologicalSort(void) {
   vector<node_t> ts;       /// Nodes in topological order
   ts.reserve(num_nodes());
   vector<node_t> d(n, 0);  /// Induced vertex incoming degree
   queue<node_t>  q;        /// Nodes to be explored yet
   /// Initialize node degrees and scan node without incoming arcs
   for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
      d[it->id] = it->in_degree();
      if ( it->in_degree() == 0 )
         q.push( it->id );
   }
   /// Core loop of topological sort
   while ( !q.empty() ) {
      node_t v = q.front();
      q.pop();
      ts.push_back(v);      
      for ( FSArcIterPair it = Nc[v].getIterFS(); it.first != it.second; ++it.first ) {
         node_t w = (it.first)->w;
         d[w]--;
         if ( d[w] == 0 )
            q.push(w);
      }
   }
   return ts;
}

/// Topological order on the reversed graph
vector<node_t>
DAG::topologicalSortBack(void) {
   vector<node_t> ts;      /// Nodes in topological order
   ts.reserve(num_nodes());
   vector<node_t> d(n, 0); /// Induced vertex incoming degree
   queue<node_t> q;        /// Nodes to be explored yet
   /// Initialize node degrees and scan node without incoming arcs
   for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
      d[it->id] = it->out_degree();
      if ( it->out_degree() == 0 )
         q.push( it->id );
   }
   /// Core loop of topological sort
   while ( !q.empty() ) {
      node_t v = q.front();
      q.pop();
      ts.push_back(v);
      for ( BSArcIterPair it = Nc[v].getIterBS(); it.first != it.second; ++it.first ) {
         node_t w = (it.first)->v;   /// It has to consider the arc as reversed
         d[w]--;
         if ( d[w] == 0 )
            q.push(w);
      }
   }
   return ts;
}

///--------------------------------------------------------------------------------
/// Solve the lagrangian subproblem with a subgradient procedure [B&C1989]
int 
DAG::subgradient(node_t Source, node_t Target, cost_t& LB, cost_t& UB) {
   /// Return value: true if it has removed at least an arc
   edge_t m0 = arcsLeft();
   /// Initialize new arc reduced cost
   ArcGradView W;
   vector<CostResources> Rf(n);
   vector<CostResources> Rb(n);
   for ( int v = 0; v < n; ++v ) {
      Rf[v].setData(0.0,k);
      Rb[v].setData(0.0,k);
   }
   
   for ( int l = 0; l < k; ++l ) 
      alpha[l] = 0.0;
   
   /// Subgradient loop
   for ( int iter = 0; iter < 10; ++iter ) {
      //for ( int l = 0; l < k; ++l ) 
         //fprintf(stdout," %.5f",alpha[l]);
      //fprintf(stdout," LAGG\n");
      /// Initialize new arc costs
      for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ++nit ) {
         for ( FSArcIterPair it = nit->getIterFS(); it.first != it.second; ++it.first ) {
            FSArcIter a = it.first;
            a->d = a->c;
            for ( int l = 0; l < k; ++l ) 
               a->d += a->r[l]*alpha[l];
         }
      }
      /// Update offset of the objective function
      UBoff = 0.0;
      for ( int l = 0; l < k; ++l ) 
         UBoff -= alpha[l]*U[l];
      /// Start double shortest path computation
      dag_ssp_all(Source, W, Rf, Pf, Df );
      /// If destination is reachable, update the lower bound
      if ( Pf[Target] == Target ) 
         return 2;

      if ( Pf[Target] != Target ) 
         LB = std::max<cost_t>(LB,Df[Target]+UBoff);
      if ( fabs(LB-ceil(LB)) > 1e-05 )
         LB = int(ceil(LB));
      else 
         LB = int(round(LB));
      /// The shortest path is not the optimum path: compute the reversed path
      dag_ssp_back_all(Target, W, Rb, Pb, Db );

      /// Filter first on the vertex and then on the forward star
      for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ) {
         node_t v = nit->id;
         if ( Df[v] + Db[v] + UBoff + EPS > UB && v != Source && v != Target) {
            NodeIter vit(nit);
            ++nit;
            clearVertex(vit);
         } else {
            filterCostFS(nit, W, Pf, Pb, Df, Db, Rf, Rb, Source, Target, UB, UBoff ); 
            ++nit; 
         }
      }
      /// Update subgradients
      cost_t Hsqr = 0.0;
      for ( int l = 0; l < k; ++l ) {
         H[l] = -U[l] + Rf[Target].r[l];
         Hsqr += H[l]*H[l];
      }
      cost_t T = f*(UB-LB)/Hsqr;

      for ( int l = 0; l < k; ++l )
         alpha[l] = std::max<cost_t>(0.0, alpha[l]+T*H[l]);
   }

   if (arcsLeft() < m0)
      return 1;
   else
      return 0;
}

///--------------------------------------------------------------------------------
void
DAG::filter( node_t Source, node_t Target, cost_t&LB, cost_t&UB ) {
   /// Iterate until an arc is removed
   subgradient(Source,Target,LB,UB);
}

/// Filter the arcs
ExecStatus DAG::filterArcs( int n, int m, ViewArray<Item>& x, Space& home) {
   for ( int i = 0; i < n; ++i ) {
      int dom_size = 0;
      int bin = -1;
      for ( IntVarValues j(x[i].bin()); j(); ++j ) {
         if ( Nc[i*m+j.val()].in_degree() > 0 ) {
            dom_size++;
            bin = j.val();
         }
      }
      if ( dom_size == 0 || (x[i].bin().assigned() && x[i].bin().val() != bin ))
         return ES_FAILED;
      if ( dom_size == 1 ) 
         GECODE_ME_CHECK(x[i].bin().eq(home, bin));
   }

   return ES_NOFIX;
}


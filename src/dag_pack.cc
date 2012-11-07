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
bool 
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
   for ( int iter = 0; iter < 20; ++iter ) {
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
      if ( Pf[Target] != Target ) 
         LB = std::max<cost_t>(LB,ceil(Df[Target]+UBoff));
      fprintf(stdout,"%.3f ",LB);
      if ( Pf[Target] == Target || LB + EPS  >= UB ) {
         LB = UB;
         return false;
      }
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
      cost_t T = f*(1.0)/Hsqr;

      for ( int l = 0; l < k; ++l )
         alpha[l] = std::max<cost_t>(0.0, alpha[l]+T*H[l]);
   }

   return (arcsLeft() < m0);
}

///--------------------------------------------------------------------------------
void
DAG::filter( node_t Source, node_t Target, cost_t&LB, cost_t&UB ) {
   bool removed = true;
   bool newUB = true;
   cost_t UB_before = UB;     

   /// Data for resource based filtering
   vector<CostResources> Rf(n);
   vector<CostResources> Rb(n);
   while ( removed || newUB ) {
      for ( int v = 0; v < n; ++v ) {
         Rf[v].setData(0.0,k);
         Rb[v].setData(0.0,k);
      }
      fprintf(stdout, "Loop %.1f %.1f\n", LB, UB);
      /// Stop conditions
      removed = false; /// MANCA IL CONTROLLO SUGLI ARCHI RIMOSSI !!!
      newUB = false;
      UB_before = UB;     

      for ( int l = 0; l < 0; ++l ) {
         fprintf(stdout, "filter resource %d\n", l);
         ArcResView W(l);
         dag_ssp_all(Source, W, Rf, Pf, Df );
         fprintf(stdout,"RES %d %f\n",U[l],Df[Target]);
         if ( Pf[Target] == Target ) {  /// If the destination is no longer reacheable
            fprintf(stdout, "[feasible] Optimal path found of value UB %.1f\n", UB);
            LB = UB;
            return;
         }
         dag_ssp_back_all(Target, W, Rb, Pb, Db );
         /// Reduce the graph
         for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ) {
            node_t v = nit->id;
            if ( v != Source && v != Target && Df[v] + Db[v] > U[l] ) {
               NodeIter vit = nit;
               ++nit;
               clearVertex(vit);
               removed = true;
            } else {
               /// Filter the arcs
               removed = filterResourceFS(nit, W, Pf, Pb, Df, Db, Rf, Rb, UB, l, Source, Target) || removed;
               ++nit;	  
            }
         }
      }
      /// Cost Based Filtering (if shortest path == UB then it is the optimum)
      if ( false ) {
         ArcCostView W;

         dag_ssp_all(Source, W,  Rf, Pf,  Df);
         /// Check for new lower bound
         if ( Pf[Target] != Target ) {
            LB = std::max(LB,computeCost(Rf[Target].c));
            /// Update the best multipliers (those that gave the lower bound)
         }
         if ( Pf[Target] == Target || LB+EPS >= UB ) {
            LB = UB;
            return;
         }
         /// Reverse shortest path tree
         dag_ssp_back_all(Target, W,  Rb, Pb,  Db);
         /// Reduce the graph
         for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ) {
            node_t v = nit->id;
            if ( v != Source && v != Target && computeCost(Df[v]+Db[v]) > UB ) {
               NodeIter vit = nit;
               ++nit;
               clearVertex(vit);
               removed = true;
            } else {
               removed = filterCostFS(nit, W, Pf, Pb, Df, Db, Rf, Rb, Source, Target, UB) || removed;
               ++nit;
            }
         }
      }
      /// Iterate until an arc is removed, or the upper bound is changed
      removed = subgradient(Source,Target,LB,UB) || removed;
      newUB = ( UB_before > UB );
   }
   fprintf(stdout,"Nodes left %d\n",(int)num_nodes());
      
}

void
DAG::printArcs(int n, int m) {
   vector<int> D(n*m+2,0);
   for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ++nit ) {
      bool flag = false;
      for ( FSArcIterPair it = nit->getIterFS(); it.first != it.second; ++it.first ) {
         fprintf(stdout,"%d %d\t", it.first->v, it.first->w);
         flag = true;
         D[it.first->w]++;
         //int i = it.first->w / m;
         //int j = it.first->w % m;
      }
      if (flag)
         fprintf(stdout,"\n");
   }
   int tot = 0;
   for ( int i = 0; i < n; ++i ) {
      int dom = 0;
      int dif = 0;
      for ( int j = 0; j < m; ++j ) {
         if ( D[i*m+j] > 0 ) {
            dom++;
            dif++;
         }
      }
      if ( dif <= 1 )
         fprintf(stdout,"#x[%d]=%d \t",i,dom);
      tot += dom;
   }
   fprintf(stdout,"\nDomTot %d\t",tot);
}

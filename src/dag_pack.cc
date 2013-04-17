#include "dag_pack.hh"
#include "cost_resources_pack.hh"
#include "path.hh"
#include "qsopt.h"

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
   
   //for ( int l = 0; l < k; ++l ) 
      //alpha[l] = 0.0;
   
   /// Subgradient loop
   for ( int iter = 0; iter < 130; ++iter ) {
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
      if ( Pf[Target] == -1 ) 
         return 2;

      if ( Pf[Target] != -1 ) 
         LB = std::max<cost_t>(LB,Df[Target]+UBoff);
      if ( fabs(LB-ceil(LB)) > EPS )
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
int
DAG::filter( node_t Source, node_t Target, cost_t& LB, cost_t& UB ) {
   /// Iterate until an arc is removed
   //return subgradient(Source, Target, LB, UB);
   return cuttingPlanes(Source, Target, LB, UB);
}

/// Filter the arcs
ExecStatus DAG::filterArcs( int n, int m, int k, ViewArray<IntView>& x, Space& home) {
   for ( int i = 0; i < n; ++i ) {
      int dom_size = 0;
      int bin = -1;
      for ( IntVarValues j(x[i]); j(); ++j ) {
         if ( Nc[i*m+j.val()].in_degree() > 0 ) {
            dom_size++;
            bin = j.val();
         }
      }
      if ( dom_size == 0 || (x[i].assigned() && x[i].val() != bin ))
         return ES_FAILED;
      if ( dom_size == 1 ) 
         GECODE_ME_CHECK(x[i].eq(home, bin));
   }

   return ES_NOFIX;
}

///------------------------------------------------------------------------------------------
/// Compute the lagrangian upper bound and the corresponding dual multipliers
int 
DAG::cuttingPlanes ( node_t Source, node_t Target, cost_t& LB, cost_t& UB ) {
   /// Return status: 0) ok 1) failed 2) filtered
   int c_status = 0;

   /// Gradient View for the separation algorithm
   ArcGradView W;
   vector<CostResources> Rf(n);
   vector<CostResources> Rb(n);
   for ( int v = 0; v < n; ++v ) {
      Rf[v].setData(0.0,k);
      Rb[v].setData(0.0,k);
   }

   /// Data structure for calling the underlying ANSI/C LP solver
   QSprob   model;
   int      error  = 0;
   int      status = 0;
   int      *cind  = NULL;
   double   *cval  = NULL;
   double   *xbar  = NULL;
   double   lb     = 0;
   int      iter   = 0;

   cind = (int*)malloc(sizeof(int) * (k+1) );
   if (!cind) goto QUIT;
   cval = (double*)malloc(sizeof(double) * (k+1) );
   if (!cval) goto QUIT;
   xbar = (double*)malloc(sizeof(double) * (k+1) );
   if (!xbar) goto QUIT;

   /// Master Problem
   model = QScreate_prob("lag_continuo", QS_MAX);
   if (model == NULL) {
      fprintf(stderr, "Error: could not create the problem\n");
      exit(EXIT_FAILURE);
   }

   /// Add the problem variables (u, v_1, ..., v_k)
   /// Add variable |u|
   error = QSnew_col(model, 1.0, -QS_MAXDOUBLE, 100*UB, (const char*) NULL);
   if (error) goto QUIT;

   /// Add variables |v|  
   for (int i = 0; i < k; ++i) {
      error = QSnew_col(model, U[i], -QS_MAXDOUBLE, 0.0, (const char*) NULL);
      if (error) goto QUIT;
   }

   /// Constraint on a path for each constraint
   for ( int l = 0; l < k ; ++l ) {
      /// Compute the shortest path using the resource consumption
      dag_ssp(Source, ArcResView(l), Pf, Df);
      if ( Pf[Target] == -1 ) {
         c_status = 2;
         goto QUIT;
      }
      /// Make the path
      Path path(*this, Source, Target, Pf);
      /// First set the variable |u|
      cind[0] = 0;
      cval[0] = 1;
      /// Then set the variables |v|
      for ( int i = 0; i < k; ++i ) {
         cind[i+1] = i+1;
         cval[i+1] = path.Rc[i];
      }
      /// Take the rhs
      double rhs = path.c;
      error = QSadd_row(model, k+1, cind, cval, rhs, 'L', (const char*) NULL);
      if (error) goto QUIT;
   }

//   error = QSwrite_prob(model, "lag.lp", "LP");
  // if (error) goto QUIT;

   QSset_param (model, QS_PARAM_DUAL_PRICING, QS_PRICE_DDANTZIG);
   do {
      /* Optimize */
      error = QSopt_dual(model, &status);
      if (error) goto QUIT;

      if (status == QS_LP_UNBOUNDED) {
         fprintf(stdout, "Unbounded\n");
         goto QUIT;
      }

      if (status == QS_LP_INFEASIBLE) {
         fprintf(stdout, "Infeasible\n");
         goto QUIT;
      }

      QSget_objval(model, &lb);
      //fprintf(stdout, "LP %f - status %d\n", lb, status);
      /// Take the current LP decision vector
      error = QSget_x_array(model, xbar);  /// xbar[0] <-> u; xbar[i+1] <-> v_i
      if (error) goto QUIT;

      double u_bar = xbar[0];

      /// Problema di separazione: Cammino Minimo su DAG
      /// Set the edge scaled weight as: c_e = c_e - v_1.r_e - ... - v_k.r_e
      for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ++nit ) {
         for ( FSArcIterPair it = nit->getIterFS(); it.first != it.second; ++it.first ) {
            FSArcIter a = it.first;
            a->d = a->c;
            for ( int l = 0; l < k; ++l ) 
               a->d -= xbar[l+1]*a->r[l];
         }
      }

      dag_ssp(Source, W, Pf, Df);
      if ( Pf[Target] == -1 ) 
         exit(EXIT_FAILURE); 

      Path path(*this, Source, Target, Pf);
      
      if ( fabs(path.d - u_bar) < EPS || path.d > u_bar + EPS ) {
         QSget_objval(model, &lb);
         goto FILTER;
      }

      /// First set the variable |u|
      cind[0] = 0;
      cval[0] = 1;
      /// Then set the variables |v|
      for ( int i = 0; i < k; ++i ) {
         cind[i+1] = i+1;
         cval[i+1] = path.Rc[i];
      }
      /// Take the rhs
      double rhs = path.c;
      error = QSadd_row(model, k+1, cind, cval, rhs, 'L', (const char*) NULL);
      if (error) goto QUIT;

   } while ( iter++ < 50 );

FILTER:
   if ( iter < 50 ) {
      /// Problema di separazione: Cammino Minimo su DAG
      /// Set the edge scaled weight as: c_e = c_e - v_1.r_e - ... - v_k.r_e
      for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ++nit ) {
         for ( FSArcIterPair it = nit->getIterFS(); it.first != it.second; ++it.first ) {
            FSArcIter a = it.first;
            a->d = a->c;
            for ( int l = 0; l < k; ++l ) 
               a->d -= xbar[l+1]*a->r[l];
         }
      }

      /// Update offset of the objective function
      UBoff = 0;
      for ( int l = 0; l < k; ++l ) 
         UBoff += xbar[l+1]*U[l];
      /// Start double shortest path computation
      dag_ssp_all(Source, W, Rf, Pf, Df );
      //fprintf(stdout,"LLBB %f \n", UBoff+Df[Target]);
      /// If destination is reachable, update the lower bound
      if ( Pf[Target] == -1 ) {
         //fprintf(stdout, "Disconnected\n");
         c_status = 2;
         goto QUIT;
      }
      //fprintf(stdout,"\n");

      LB = std::max<cost_t>(LB,Df[Target]+UBoff);

      //if ( fabs(LB-ceil(LB)) > EPS )
         //LB = int(ceil(LB));
      //else 
         //LB = int(round(LB));
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
   }
   else {
      fprintf(stdout,"WARNING\n");
      c_status = 1;
      //exit(0);
   }
   
QUIT:
   free(cind);
   free(cval);
   free(xbar);

   QSfree_prob(model);

   if ( c_status == 1 )
      LB = std::max<cost_t>(lb, LB);

   return c_status;
}

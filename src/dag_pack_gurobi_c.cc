#include "dag_pack.hh"
#include "cost_resources_pack.hh"
#include "path.hh"

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
   /// Initialize new arc reduced cost
   ArcGradView W;
   for ( int v = 0; v < n; ++v ) {
      Rf[v].setData(0.0,k);
      Rb[v].setData(0.0,k);
   }
   
   //for ( int l = 0; l < k; ++l ) 
      //alpha[l] = 0.0;
   
   /// Subgradient loop
   int iter;
   for ( iter = 0; iter < 50; ++iter ) {
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
         if ( Df[v] + Db[v] + UBoff - EPS > UB && v != Source && v != Target) {
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

   if ( iter >= 49 )
      return 1;
   return 0;
}

///--------------------------------------------------------------------------------
int
DAG::filter( node_t Source, node_t Target, cost_t& LB, cost_t& UB ) {
   /// Iterate until an arc is removed
//   if ( pool.empty() )
      return cuttingPlanes(Source, Target, LB, UB);
  // else
    //  return subgradient(Source, Target, LB, UB);
}

/// Filter the arcs
ExecStatus DAG::filterArcs( int n, int m, int k, ViewArray<IntView>& x, Space& home) {
   for ( int i = 0; i < n; ++i ) {
      //fprintf(stdout,"\n");
      if ( !x[i].assigned() ) {
         int dom_size = 0;
         int bom_size = 0;
         int bin = -1;
         for ( IntVarValues j(x[i]); j(); ++j ) {
            //fprintf(stdout,"x[%d]=%d\t", i, j.val());
            if ( Nc[i*m+j.val()].in_degree() > 0 ) {
               dom_size++;
               bin = j.val();
            }
            if ( Nc[i*m+j.val()].out_degree() > 0 ) {
               bom_size++;
               bin = j.val();
            }
         }
         if ( dom_size == 0 || bom_size == 0) { // || (x[i].assigned() && x[i].val() != bin ))
            return ES_FAILED;
         }

         if ( dom_size == 1 || bom_size == 1 ) 
            GECODE_ME_CHECK(x[i].eq(home, bin));
      }
   }

   return ES_NOFIX;
}

///------------------------------------------------------------------------------------------
/// Compute the lagrangian upper bound and the corresponding dual multipliers
int 
DAG::cuttingPlanes ( node_t Source, node_t Target, cost_t& LB, cost_t& UB ) {
   /// Return status: 0) ok 1) iter limit 2) failed
   int c_status = 0;

   /// Gradient View for the separation algorithm
   ArcGradView W;
   for ( int v = 0; v < n; ++v ) {
      Rf[v].setData(0.0,k);
      Rb[v].setData(0.0,k);
   }

   /// Data structure for calling the underlying ANSI/C LP solver
   int      error  = 0;
   int      status = 0;
   double   lb     = 0;
   int      iter   = 0;
   double   t0 = 0;
   double   t1 = 0;
   double   t_lp = 0;
   int      nz = 0;

   timer TT;

   if ( pool.empty() ) {
      /// Constraint on a path for each constraint
      for ( int l = 0; l < k ; ++l ) {
         /// Compute the shortest path using the resource consumption
         dag_ssp(Source, ArcResView(l), Pf, Df);
         if ( Pf[Target] == -1 ) {
            c_status = 2;
            goto QUIT;
         }
         /// Make the path
         Path* path = new Path(*this, Source, Target, Pf);
         /// First set the variable |u|
         cind[0] = 0;
         cval[0] = 1;
         nz = 1;
         /// Then set the variables |v|
         for ( int i = 0; i < k; ++i ) {
            if ( path->Rc[i] > 0 ) {
               cind[nz] = i+1;
               cval[nz] = path->Rc[i];
               nz++;
            }
         }
         /// Take the rhs
         double rhs = path->c;
         error = GRBaddconstr(model, nz, cind, cval, GRB_LESS_EQUAL, rhs, (const char*) NULL);
         if (error) goto QUIT;
         pool.push_back(path);
      }
   } else {
      int np = pool.size();
      for ( int h = 0; h < np; ++h ) {
         /// Take the rhs
         double rhs = pool[h].updateCost(*this);
         error = GRBsetdblattrelement(model, GRB_DBL_ATTR_RHS, h, rhs);
         if (error) goto QUIT;
      }
   } 

   t0 = TT.elapsed();

   do {
      /* Optimize */
      t1 = TT.elapsed();
      error = GRBoptimize(model);//QSopt_dual(model, &status);
      t_lp += TT.elapsed()-t1;
      if (error) goto QUIT;

      error = GRBgetintattr(model, "Status", &status);
      if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE ) {
         fprintf(stdout, "Unbounded or Infeasible\n");
         c_status = 1;
         goto QUIT;
      }

      /// Take the current LP decision vector
      error = GRBgetdblattrarray(model, "X", 0, k+1, xbar);
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

      //dag_ssp(Source, W, Pf, Df);

      //if ( Pf[Target] == -1 ) {
         //c_status = 2;
         //goto QUIT;
      //}

      /// Start double shortest path computation
      dag_ssp_all(Source, W, Rf, Pf, Df );
      
      /// If destination is reachable, update the lower bound
      if ( Pf[Target] == -1 ) {
         c_status = 2;
         goto QUIT;
      }
      
      /// Update offset of the objective function
      UBoff = 0;
      for ( int l = 0; l < k; ++l ) 
         UBoff += xbar[l+1]*U[l];
      LB = std::max<cost_t>(LB,Df[Target]+UBoff);
      if ( fabs(LB-ceil(LB)) > EPS )
         LB = int(ceil(LB));
      else 
         LB = int(round(LB));

      if ( LB > UB ) {
         //fprintf(stdout, "Disconnected\n");
         c_status = 2;
         goto QUIT;
      }

      /// The shortest path is not the optimum path: compute the reversed path
      dag_ssp_back_all(Target, W, Rb, Pb, Db );

      /// Filter first on the vertex and then on the forward star
      for ( NodeIter nit = N.begin(), nit_end = N.end(); nit != nit_end; ) {
         node_t v = nit->id;
         if ( Df[v] + Db[v] + UBoff - EPS > UB && v != Source && v != Target) {
            NodeIter vit(nit);
            ++nit;
            clearVertex(vit);
         } else {
            filterCostFS(nit, W, Pf, Pb, Df, Db, Rf, Rb, Source, Target, UB, UBoff ); 
            ++nit; 
         }
      }

      Path* path = new Path(*this, Source, Target, Pf);
      
      if ( fabs(path->d - u_bar) < EPS || path->d > u_bar + EPS ) {
         GRBgetdblattr(model, "ObjVal", &lb);
         goto QUIT;
      }

      /// First set the variable |u|
      cind[0] = 0;
      cval[0] = 1;
      nz = 1;
      /// Then set the variables |v|
      for ( int i = 0; i < k; ++i ) {
         if ( path->Rc[i] > 0 ) {
            cind[nz] = i+1;
            cval[nz] = path->Rc[i];
            nz++;
         }
      }
      /// Take the rhs
      double rhs = path->c;
      error = GRBaddconstr(model, nz, cind, cval, GRB_LESS_EQUAL, rhs, (const char*) NULL);
      if (error) goto QUIT;
      pool.push_back(path);
   } while ( iter++ < 200 );

   if ( iter >= 200 ) 
      //fprintf(stdout,"WARNING\n");
      c_status = 1;
   
QUIT:
   if ( c_status == 0 )
      LB = std::max<cost_t>(lb, LB);
   for ( int l = 0; l < k; ++l )
      alpha[l] = fabs(xbar[l+1]);
   //fprintf(stdout, "TiFilter %.5f %.5f %.5f %d\n", TT.elapsed(), TT.elapsed()-t0, t_lp, iter); 
   return c_status;
}

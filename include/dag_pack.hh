#ifndef __MY_DAG_PACK
#define __MY_DAG_PACK

/// From STL library
#include <vector>
using std::vector;

#include <stack>
using std::stack;

#include <queue>
using std::queue;

/// Numerics limits
#include <limits>
using std::numeric_limits;

/// My includes
#include "arc.hh"
#include "node.hh"
#include "cost_resources_pack.hh"

#include <stdio.h>

#include <gecode/int.hh>
using namespace Gecode;
using namespace Gecode::Int;

#include "cost_binpacking.hh"
using namespace Gecode::Int::CostBinPacking;

/// Constants
const cost_t EPS = 1e-09;

///--------------------------------------------------------------------------------
/// Class of graph to compute RCSP with superadditive cost
class DAG {
   private:
      const node_t  n;
      edge_t  m;
      const int     k;

      vector<Arc>     A;    /// Arcs  container 
      vector<Node>    Nc;   /// Nodes container

      vector<node_t>  Pf;   /// Forward  predecessor vector
      vector<node_t>  Pb;   /// Backward predecessor vector
      vector<dist_t>  Df;   /// Forward  distance vector
      vector<dist_t>  Db;   /// Backward distance vector

      /// Initialize distance vector with Infinity
      /// Maybe it is better to intialize with an upper bound on the optimal path (optimal rcsp path)
      const dist_t Inf;
      
      cost_t UBoff;
      vector<cost_t> alpha;
      vector<cost_t> H;
      cost_t H0;
      cost_t zb;
      cost_t f;

   public:
      resources     U;  /// Upper limits for the *k* resources
      NodeList      N;  /// List view (Intrusive list)

      ///Standard constructor
      DAG ( node_t _n, edge_t _m, const resources& _U ) 
         : n(_n), m(_m), k(_U.size()), Nc(n), Pf(n), Pb(n), Df(n), Db(n),
         Inf(std::numeric_limits<dist_t>::max()),
         alpha(k,0.0), H(k,0.0), U(_U)
   {
      assert( n < std::numeric_limits<node_t>::max()  &&
            m < std::numeric_limits<edge_t>::max() );
      /// Reserve memory
      A.reserve(m);
      m = 0;
      /// Initialize the node set
      for ( node_t i = 0; i < n; ++i ) {
         Nc[i].setData(i);
         N.push_back( Nc[i] );
      }
      /// Subgradient
      H0 = 0.0;
      zb = 0.0;
      f  = 2;
   }

      /// Basic getters
      inline node_t num_nodes(void) const { return N.size(); }
      inline edge_t num_arcs(void)  const { return A.size(); }

      /// Basic setter: called only by the constructor
      void addArc( node_t i, node_t j, cost_t c, const resources& r ) {    
         edge_t arc_id = static_cast<edge_t>(A.size());
         /// Added it once for ever
         A.push_back( Arc() );
         A[arc_id].setData( arc_id, i, j, c, r );
         m++;
         /// Allocated space once for ever and use view to iterate on them
         Nc[i].addForwArc(A[arc_id]);
         Nc[j].addBackArc(A[arc_id]);
      }

      inline void removeArc( FSArcIter a ) {
         //fprintf(stdout,"removed %d %d\n", a->v, a->w);
         Nc[a->v].removeForwArc( a );
         Arc& b = A[a->id];
         Nc[a->w].removeBackArc( b );
      }

      inline void removeOppositeArc( FSArcIter a ) {
         Arc& b = A[a->id];
         Nc[a->w].removeBackArc( b );
      }

      inline void removeOppositeArc( BSArcIter a ) {
         Arc& b = A[a->id];
         Nc[a->v].removeForwArc( b );
      }

      edge_t arcsLeft(void) {
         edge_t left = 0;
         for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) 
            left += it->out_degree();
         return (m-left);
      }
      
      ///--------------------------------------------------------------------------------
      /// Test if the given graph is acyclic
      bool isAcyclic (void) {
         return (static_cast<int>(topologicalSort().size()) == n);
      }

      ///--------------------------------------------------------------------------------
      /// Shortest Path for a DAG: record only predecessor and length of each node
      template <class LengthMap>
         void
         dag_ssp ( node_t S, const LengthMap& W,
               vector<node_t>& P, vector<dist_t>& D ) {    
            /// Initialize predecessor and distance vectors
            for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
               P[it->id] = it->id;
               D[it->id] = Inf;
            }
            /// Initialize the source distance
            D[S] = 0;

            /// Get the topological sort 
            vector<node_t> ts = topologicalSort();

            /// Start the algorithm (*Reaching*, see Denardo and Fox 1979, JSTOR)
            vector<node_t>::const_iterator it = ts.begin();
            vector<node_t>::const_iterator it_end = ts.end();
            while ( it != it_end ) {
               node_t u = *it++;
               if ( D[u] < Inf ) { /// If D[u]=Inf then D[u]+Wuv is Inf as well
                  for ( FSArcIterPair it = Nc[u].getIterFS(); it.first != it.second; ++it.first ) {
                     node_t v = (it.first)->w;
                     dist_t Wuv = static_cast<dist_t>(W(*it.first));
                     if ( D[v] > D[u] + Wuv ) {
                        D[v] = D[u] + Wuv; 
                        P[v] = u;
                     }
                  }
               }
            }
         }

      ///--------------------------------------------------------------------------------
      /// Shortest Path for a DAG
      template <class LengthMap>
         void
         dag_ssp_back ( node_t S, const LengthMap& W,
               vector<node_t>& P, vector<dist_t>& D ) {
            /// Initialize predecessor vector
            for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
               P[it->id] = it->id;
               D[it->id] = Inf;
            }
            /// Init the source node
            D[S] = 0;
            /// Get the topological sort 
            vector<node_t> ts = topologicalSortBack();
            /// Start the algorithm
            vector<node_t>::const_iterator it = ts.begin();
            vector<node_t>::const_iterator it_end = ts.end();
            while ( it != it_end ) {
               node_t u = *it++;
               if ( D[u] < Inf ) { /// If D[u]=Inf then D[u]+Wuv is Inf as well
                  for ( BSArcIterPair it = Nc[u].getIterBS(); it.first != it.second; ++it.first ) {
                     node_t v = (it.first)->v;   /// It has to consider the arc as reversed
                     dist_t Wuv = static_cast<dist_t>(W(*it.first));
                     if ( D[v] > D[u] + Wuv ) {
                        D[v] = D[u] + Wuv; /// Arc Length
                        P[v] = u;
                     }
                  }
               }
            }
         }

      ///--------------------------------------------------------------------------------
      /// Shortest Path for a DAG: record only predecessor and length of each node
      template <class LengthMap>
         void
         dag_ssp_all ( node_t S, const LengthMap& W, vector<CostResources>& R,
               vector<node_t>& P, vector<dist_t>& D ) {
            /// Initialize predecessor vector
            for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
               P[it->id] = it->id;
               D[it->id] = Inf;
            }
            /// Init the source node
            D[S] = 0;
            /// Get the topological sort 
            vector<node_t> ts = topologicalSort();    
            /// Start the algorithm (*Reaching*, see Denardo and Fox 1979, JSTOR)
            vector<node_t>::const_iterator it = ts.begin();
            vector<node_t>::const_iterator it_end = ts.end();
            while ( it != it_end ) {
               node_t u = *it++;
               if ( D[u] < Inf ) { /// If D[u]=Inf then D[u]+Wuv is Inf as well
                  for ( FSArcIterPair it = Nc[u].getIterFS(); it.first != it.second; ++it.first ) {
                     node_t v = (it.first)->w; /// Arc length (it depends on the LengthView)
                     dist_t Wuv = static_cast<dist_t>(W(*it.first));
                     if ( D[v] > D[u] + Wuv ) {
                        D[v] = D[u] + Wuv; 
                        P[v] = u;
                        R[v].add(R[u], *it.first);
                     }
                  }
               }
            }    
         }

      ///--------------------------------------------------------------------------------
      /// Shortest Path for a reversed DAG: record only predecessor and length of each node
      template <class LengthMap>
         void
         dag_ssp_back_all ( node_t S, const LengthMap& W, vector<CostResources>& R,
               vector<node_t>& P, vector<dist_t>& D ) {
            /// Initialize predecessor vector
            for ( NodeIter it = N.begin(), it_end = N.end(); it != it_end; ++it ) {
               P[it->id] = it->id;
               D[it->id] = Inf;
            }
            /// Init the source node
            D[S] = 0;
            /// Get the topological sort 
            vector<node_t> ts = topologicalSortBack();
            /// Start the algorithm
            vector<node_t>::const_iterator it = ts.begin();
            vector<node_t>::const_iterator it_end = ts.end();
            while ( it != it_end ) {
               node_t u = *it++;
               if ( D[u] < Inf ) { /// If D[u]=Inf then D[u]+Wuv is Inf as well
                  for ( BSArcIterPair it = Nc[u].getIterBS(); it.first != it.second; ++it.first ) {
                     node_t v = (it.first)->v;   /// It has to consider the arc as reversed
                     dist_t Wuv = static_cast<dist_t>(W(*it.first));
                     if ( D[v] > D[u] + Wuv ) {
                        D[v] = D[u] + Wuv; /// Arc Length
                        P[v] = u;
                        R[v].add(R[u], *it.first);
                     }
                  }
               }
            }
         }

      /// Filter the forward star of a vertex
      template <class LengthMap>
         bool filterResourceFS(NodeIter vit, const LengthMap& W,
               const vector<node_t>& Pf, const vector<node_t>& Pb,
               const vector<dist_t>& Df, const vector<dist_t>& Db,
               const vector<CostResources>& Rf, const vector<CostResources>& Rb,
               cost_t& UB, int l, node_t Source, node_t Target) 
         {
            bool removed = false;
            node_t v = vit->id;
            FSArcIterPair it = vit->getIterFS();
            /// Filter the arcs
            while ( it.first != it.second ) {
               FSArcIter a = it.first;
               ++it.first;  /// First increment the pointer, for safe in-place removals
               node_t w = a->w; /// Arc length (it depends on the LengthView)
               dist_t Wuv = static_cast<dist_t>(W(*a));
               if ( (v != Source && Pf[v] == v) || (w != Target && Pb[w] == w) || 
                     Df[v] + Wuv + Db[w] > U[l] ) {
                  removeArc(a);
                  removed = true;
               }
            }
            return removed;
         }

      /// Filter the forward star of a vertex with the cost
      template <class LengthMap>
         bool filterCostFS(NodeIter vit, const LengthMap& W,
               const vector<node_t>& Pf, const vector<node_t>& Pb,
               const vector<dist_t>& Df, const vector<dist_t>& Db,
               const vector<CostResources>& Rf, const vector<CostResources>& Rb,
               node_t Source, node_t Target, cost_t& UB, cost_t UB_off=EPS ) 
         {
            bool removed = false;
            node_t v = vit->id;
            FSArcIterPair it = vit->getIterFS();
            /// Filter the arcs
            while ( it.first != it.second ) {
               FSArcIter a = it.first;
               ++it.first;  /// First increment the pointer, for safe in-place removals
               node_t w = a->w; /// Arc length (it depends on the LengthView)
               dist_t Wuv = static_cast<dist_t>(W(*a));
               if ( (v != Source && Pf[v] == v) || (w != Target && Pb[w] == w) || 
                     computeCost(Df[v]+Wuv+Db[w]+UB_off) > UB ) {
                  removeArc(a);
                  removed = true;
               }
            }
            return removed;
         }

      /// Methods implemented in the .cc file
      vector<node_t> topologicalSort(void);
      vector<node_t> topologicalSortBack(void);

      void clearVertex ( NodeIter vit );
      bool subgradient ( node_t Source, node_t Target, cost_t& LB, cost_t& UB );
      void filter      ( node_t Source, node_t Target, cost_t& LB, cost_t& UB );
      void printArcs   ( int n, int m ); 
      void filterArcs  ( int n, int m, ViewArray<Item>& x);
}; /// End intrusive graph class

#endif /// __MY_DAG_PACK

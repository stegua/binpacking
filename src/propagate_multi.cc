#include "cost_multibin.hh"
#include <stdio.h>
#include <dag_pack.hh>
#include <map>
using std::map;

namespace Gecode { namespace Int { namespace CostMultiBinPacking {
   PropCost MultiPack::cost(const Space&, const ModEventDelta&) const {
      return PropCost::cubic(PropCost::HI,x.size());
   }

   Actor* MultiPack::copy(Space& home, bool share) {
      return new (home) MultiPack(home,share,*this);
   }

  ExecStatus 
  MultiPack::propagate(Space& home, const ModEventDelta& med) {
    /// Check if SUBSUMED
    {
       int i = 0;
       while ( i<n && x[i].assigned()) { i++; };
       if ( i == x.size() )
          return home.ES_SUBSUMED(*this);
    }

    /// Create the graph and propagates: x = x
    int N = n*m+2;
    int M = (n-1)*(m*m)+2*m;
    int K = k*m;
    int S = N-2;
    int T = N-1;

    resources U(K,0);
    for ( int l = 0; l < K; ++l )
       U[l] = y[l].max();

    DAG G (N, M, U);
    
    /// check the ammisibility
    resources B(K,0); 
    for ( int i = 0; i < n; ++i ) {
       if ( x[i].assigned() ) {
          int j = x[i].val();
          for ( int l = 0; l < k; ++l ) 
             B[j+l*m] += D[i*k+l];
       }
    }
    for ( int l = 0; l < K; ++l ) {
       if ( B[l] > y[l].max() ) {
          fprintf(stdout,"Precheck\n");
          return ES_FAILED;
       }
       else 
          GECODE_ME_CHECK(y[l].gq(home, B[l]));
    }
    
    /// Build the Directed Acyclic Graph
    typedef std::pair<int, int> NodePair;
    map< NodePair, edge_t > As;
    for ( int i = 0; i < n-1; ++i ) {
       for ( IntVarValues j(x[i+1]); j(); ++j ) {
          resources R(K,0);
          for ( int l = 0; l < k; ++l )
             R[j.val()+l*m] = D[(i+1)*k+l];//A[i+1][l];
          for ( IntVarValues h(x[i]); h(); ++h ) {
             edge_t arc_id = G.addArc( i*m+h.val(), (i+1)*m+j.val(), 0, R );
             As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))] = arc_id;
          }
       }
    }
    /// Arcs from the source node
    for ( IntVarValues j(x[0]); j(); ++j ) {
       resources R(K,0);
       for ( int l = 0; l < k; ++l )
          R[j.val()+l*m] = D[l];//x[0].size();
       edge_t arc_id = G.addArc( S, j.val(), 0, R );
       As[make_pair(S, j.val())] = arc_id;
    }
    /// Arcs to the destination node
    for ( IntVarValues j(x[n-1]); j(); ++j ) {
       resources R(K,0);
       G.addArc( (n-1)*m+j.val(), T, 0, R );
    }

    /// Filter the arcs
    cost_t LB;
    cost_t UB;
    int status = 0;
    for ( int q = 0; q < m; ++q ) {
       for ( int l = 0; l < k; ++l ) {
         LB = 0;
         UB = y[l*m+q].max();
         //fprintf(stdout,"%.1f ", UB);
         /// Archi da item a item
         for ( int i = 0; i < n-1; ++i ) {
            for ( IntVarValues j(x[i]); j(); ++j ) {
               for ( IntVarValues h(x[i+1]); h(); ++h ) {
                  if ( j.val() == q )
                     G.setArcCost(As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))], D[(i+1)*k+l]);
                  else
                     G.setArcCost(As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))], 0);
               }
            }
         }
         /// Archi dalla sorgente 
         for ( IntVarValues j(x[0]); j(); ++j ) {
            if ( j.val() == q )
               G.setArcCost( As[make_pair(S, j.val())], D[l]);
            else
               G.setArcCost( As[make_pair(S, j.val())], 0);
         }

         status = G.filter(S,T,LB,UB);

         if ( status == 2 ) {
            return ES_FAILED;
         }

         /// Aumenta il lower bound della variable di load
         if ( status == 0 ) {
            int LB0 = int(ceil(LB-0.5));
            fprintf(stdout,"#%d ", LB0);
            if ( y[l*m+q].max() < LB0 ) {
               fprintf(stdout,"\n");
               return ES_FAILED;
            }
            if ( y[l*m+q].min() < LB0 ) 
               GECODE_ME_CHECK(y[l*m+q].gq(home,LB0));
         }
       }
       fprintf(stdout,"\n");
    }
    
    if ( ES_FAILED == G.filterArcs(n,m,k,x,home) )
       return ES_FAILED;
    
    return ES_FIX;
  }

  ExecStatus
  MultiPack::post(Home home, int n, int m, int k,
      ViewArray<IntView>& y, ViewArray<IntView>& x, const IntSharedArray& D) {
    if (x.size() == 0) {
      // No items to be packed
      for (int i=y.size(); i--; )
        GECODE_ME_CHECK(y[i].eq(home,0));
      return ES_OK;
    } else if (y.size() == 0) {
      // No bins available
      return ES_FAILED;
    } else {
      // Constrain bins 
      for (int i=x.size(); i--; ) {
        GECODE_ME_CHECK(x[i].gq(home,0));
        GECODE_ME_CHECK(x[i].le(home,y.size()));
      }
      (void) new (home) MultiPack(home,n,m,k,y,x,D);
      return ES_OK;
    }
    /// Check also the z variable!
  }
}}}

// STATISTICS: int-prop


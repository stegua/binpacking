#include "cost_binpacking.hh"
#include <stdio.h>
#include <dag_pack.hh>

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
       while ( i<n && x[i].bin().assigned()) { i++; };
       if ( i == x.size() )
          return home.ES_SUBSUMED(*this);
    }

    /// Create the graph and propagates: x = x
    resources U(m, 0);
    for ( int i = 0; i < m; ++i )
       U[i] = resource_t(y[i].max());

    node_t N = n*m+2;
    edge_t M = n*m*m;
    node_t S = n*m;
    node_t T = n*m+1;
    
    DAG G (N, M, U);

    /// Build the Directed Acyclic Graph
    for ( int i = 0; i < n-1; ++i ) {
       for ( IntVarValues j(x[i].bin()); j(); ++j ) {
          for ( IntVarValues h(x[i+1].bin()); h(); ++h ) {
             if ( !(j.val() == h.val() && x[i].size() + x[i+1].size() > y[h.val()].max() ) ) {
                cost_t c = D[(i+1)*m+h.val()];
                resources R(m,0);
                R[h.val()] = x[i+1].size();
                G.addArc( i*m+j.val(), (i+1)*m+h.val(), c, R );
             }
          }
       }
    }
    /// Arcs from the source node
    for ( IntVarValues h(x[0].bin()); h(); ++h ) {
       if ( x[0].size() <= l[h.val()].max() ) {
          cost_t c = D[h.val()];
          resources R(m,0);
          R[h.val()] = x[0].size();
          G.addArc( S, h.val(), c, R );
       }
    }
    for ( IntVarValues h(x[n-1].bin()); h(); ++h ) {
       if ( x[n-1].size() <= y[h.val()].max() ) {
          cost_t c = 0.0;
          resources R(m,0);
          G.addArc( (n-1)*m+h.val(), T, c, R );
       }
    }

    /// Filter the arcs
    int status = G.subgradient(S,T,LB,UB);
   
    if ( status == 2 )
       return ES_FAILED;

    if ( z.min() < int(LB) ) 
       GECODE_ME_CHECK(z.gq(home,int(LB)));

    if ( ES_FAILED == G.filterArcs(n,m,x,home) )
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
        GECODE_ME_CHECK(x[i].bin().gq(home,0));
        GECODE_ME_CHECK(x[i].bin().le(home,y.size()));
      }
      (void) new (home) MultiPack(home,n,m,k,y,x,D);
      return ES_OK;
    }
    /// Check also the z variable!
  }
}}}

// STATISTICS: int-prop


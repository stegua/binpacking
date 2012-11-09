#include "cost_binpacking.hh"
#include <stdio.h>
#include <dag_pack.hh>

namespace Gecode { namespace Int { namespace CostBinPacking {

   /*
    *    * Packing propagator
    *       *
    *          */

   PropCost
      Pack::cost(const Space&, const ModEventDelta&) const {
         return PropCost::crazy(PropCost::HI,x.size());
      }

   Actor*
      Pack::copy(Space& home, bool share) {
         return new (home) Pack(home,share,*this);
      }

   /// Record tell information
   class TellCache {
      protected:
         /// Values (sorted) to be pruned from view
         int* _nq;
         /// Number of values to be pruned
         int _n_nq;
         /// Value to which view should be assigned
         int _eq;
      public:
         /// Initialize cache for at most \a m values
         TellCache(Region& region, int m);
         /// Record that view must be different from \a j
         void nq(int j);
         /// Record that view must be equal to \a j, return false if not possible
         void eq(int j);
         /// Perform tell to view \a x and reset cache
         ExecStatus tell(Space& home, IntView x);
   };


  forceinline
  TellCache::TellCache(Region& region, int m) 
    : _nq(region.alloc<int>(m)), _n_nq(0), _eq(-1) {}
  forceinline void 
  TellCache::nq(int j) {
    _nq[_n_nq++] = j;
  }
  forceinline void
  TellCache::eq(int j) {
    // For eq: -1 mean not yet assigned, -2 means failure, positive means value
    if (_eq == -1)
      _eq = j;
    else
      _eq = -2;
  }
  ExecStatus
  TellCache::tell(Space& home, IntView x) {
    if (_eq == -2) {
      return ES_FAILED;
    } else if (_eq >= 0) {
      GECODE_ME_CHECK(x.eq(home,_eq));
    }
    Iter::Values::Array nqi(_nq, _n_nq);
    GECODE_ME_CHECK(x.minus_v(home, nqi));
    _n_nq=0; _eq=-1;
    return ES_OK;
  }


  /*
   * Propagation proper
   *
   */
  ExecStatus 
  Pack::propagate(Space& home, const ModEventDelta& med) {
    // Number of items
    int n = x.size();
    // Number of bins
    int m = l.size();

    /// Create the graph and propagates: x = x
    cost_t LB = cost_t(z.min());
    cost_t UB = cost_t(z.max());

    resources U(m, 0);
    for ( int i = 0; i < m; ++i )
       U[i] = resource_t(l[i].max());

    node_t N = n*m+2;
    edge_t M = n*m*m;
    node_t S = n*m;
    node_t T = n*m+1;
    
    DAG G (N, M, U);

    /// Build the Directed Acyclic Graph
    for ( int i = 0; i < n-1; ++i ) 
       for ( IntVarValues j(x[i].bin()); j(); ++j ) {
          for ( IntVarValues h(x[i+1].bin()); h(); ++h ) {
             if ( !(j.val() == h.val() && x[i].size() + x[i+1].size() > l[h.val()].max() ) ) {
                cost_t c = D[(i+1)*m+h.val()];
                resources R(m,0);
                R[h.val()] = x[i+1].size();
                G.addArc( i*m+j.val(), (i+1)*m+h.val(), c, R );
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
       if ( x[n-1].size() <= l[h.val()].max() ) {
          cost_t c = 0.0;
          resources R(m,0);
          G.addArc( (n-1)*m+h.val(), T, c, R );
       }
    }
    /// Filter the arcs
    G.filter(S,T,LB,UB);
    if ( fabs(LB-ceil(LB)) > 1e-05 )
       LB = int(ceil(LB));
    else 
       LB = int(round(LB));
    
    if ( z.min() < int(LB) ) {
       //fprintf(stdout, "_ %d %.2f %.2f\n", z.min(), LB, UB);
       GECODE_ME_CHECK(z.gq(home,int(LB)));
    }

    if ( G.filterArcs(n,m,x,home) == ES_FAILED )
       return ES_FAILED;
   
    if ( LB >= UB )
       return home.ES_SUBSUMED(*this);

    return ES_FIX;
  }

  ExecStatus
  Pack::post(Home home, ViewArray<OffsetView>& l, ViewArray<Item>& x, IntView& z, const IntSharedArray& D) {
    if (x.size() == 0) {
      // No items to be packed
      for (int i=l.size(); i--; )
        GECODE_ME_CHECK(l[i].eq(home,0));
      return ES_OK;
    } else if (l.size() == 0) {
      // No bins available
      return ES_FAILED;
    } else {
      int s = 0;
      // Constrain bins 
      for (int i=x.size(); i--; ) {
        s += x[i].size();
        GECODE_ME_CHECK(x[i].bin().gq(home,0));
        GECODE_ME_CHECK(x[i].bin().le(home,l.size()));
      }
      // Constrain load
      for (int j=l.size(); j--; ) {
        GECODE_ME_CHECK(l[j].gq(home,0));
        GECODE_ME_CHECK(l[j].lq(home,s));
      }
      (void) new (home) Pack(home,l,x,z,D);
      return ES_OK;
    }
    /// Check also the z variable!
  }

}}}

// STATISTICS: int-prop


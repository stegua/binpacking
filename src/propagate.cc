#include "cost_binpacking.hh"
#include <stdio.h>


namespace Gecode { namespace Int { namespace CostBinPacking {

   /*
    *    * Packing propagator
    *       *
    *          */

   PropCost
      Pack::cost(const Space&, const ModEventDelta&) const {
         return PropCost::crazy(PropCost::HI,bs.size());
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
    int n = bs.size();
    // Number of bins
    int m = l.size();

    /// Create the graph and propagates: bs = x
    
    return ES_NOFIX;
  }

  ExecStatus
  Pack::post(Home home, ViewArray<OffsetView>& l, ViewArray<Item>& bs) {
    if (bs.size() == 0) {
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
      for (int i=bs.size(); i--; ) {
        s += bs[i].size();
        GECODE_ME_CHECK(bs[i].bin().gq(home,0));
        GECODE_ME_CHECK(bs[i].bin().le(home,l.size()));
      }
      // Constrain load
      for (int j=l.size(); j--; ) {
        GECODE_ME_CHECK(l[j].gq(home,0));
        GECODE_ME_CHECK(l[j].lq(home,s));
      }
      (void) new (home) Pack(home,l,bs);
      return ES_OK;
    }
  }

}}}

// STATISTICS: int-prop


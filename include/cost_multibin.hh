#ifndef __MY_PACK_PROPGATOR
#define __MY_PACK_PROPGATOR

#include <gecode/int.hh>


namespace Gecode { namespace Int { namespace CostMultiBinPacking {
   class MultiPack : public Propagator {
      protected:
         int n;
         int m;
         int k;
         /// Views for load of bins
         ViewArray<IntView> y;
         /// Items with bin and size
         ViewArray<IntView> x;
         /// Cost of assign processi i to machine j
         IntSharedArray  D;
         /// Constructor for posting
         MultiPack(Home home, int n, int m, int k, ViewArray<IntView>& y, ViewArray<IntView>& x, const IntSharedArray& D);
         /// Constructor for cloning \a p
         MultiPack(Space& home, bool share, MultiPack& p);
      public:
         /// Post propagator for loads \a l and items \a x
         GECODE_INT_EXPORT 
            static ExecStatus post(Home home, 
                  int n, int m, int k,
                  ViewArray<IntView>&    y, 
                  ViewArray<IntView>&    x, 
                  const IntSharedArray&  D);
         /// Perform propagation
         GECODE_INT_EXPORT 
            virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
         /// Cost function
         GECODE_INT_EXPORT 
            virtual PropCost cost(const Space& home, const ModEventDelta& med) const;
         /// Copy propagator during cloning
         GECODE_INT_EXPORT 
            virtual Actor* copy(Space& home, bool share);
         /// Destructor
         virtual size_t dispose(Space& home);
   };
}}

GECODE_INT_EXPORT void
cost_multibin(Home home, 
      int n, int m, int k,
      const IntVarArgs&      y, 
      const IntVarArgs&      x, 
      const IntSharedArray&  D,
      IntConLevel icl=ICL_DEF);
}

#include "propagate_multibin.hh"

#endif

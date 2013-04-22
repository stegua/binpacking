namespace Gecode { namespace Int { namespace CostMultiBinPacking {
   forceinline
      MultiPack::MultiPack(Home home, int n0, int m0, int k0,
            ViewArray<IntView>& y0, ViewArray<IntView>& x0, const IntSharedArray& D0, const IntSharedArray& B0 )
      : Propagator(home), n(n0), m(m0), k(k0), y(y0), x(x0), D(D0), B(B0) {
         x.subscribe(home,*this,PC_INT_DOM);
         home.notice(*this,AP_DISPOSE);
      }

   forceinline
      MultiPack::MultiPack(Space& home, bool shared, MultiPack& p) 
      : Propagator(home,shared, p), n(p.n), m(p.m), k(p.k) {
         y.update(home, shared, p.y);
         x.update(home, shared, p.x);
         D.update(home, shared, p.D);
         B.update(home, shared, p.B);
      }

   forceinline size_t 
      MultiPack::dispose(Space& home) {
         x.cancel(home,*this,PC_INT_DOM);
         home.ignore(*this,AP_DISPOSE);
         D.~IntSharedArray();
         B.~IntSharedArray();
         (void) Propagator::dispose(home);
         return sizeof(*this);
      }
}}}

// STATISTICS: int-prop


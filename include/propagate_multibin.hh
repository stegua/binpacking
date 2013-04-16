namespace Gecode { namespace Int { namespace CostMultiBinPacking {
   forceinline
      MultiPack::MultiPack(Home home, int n, int m, int k,
            ViewArray<IntView>& y0, ViewArray<IntView>& x0, const IntSharedArray& D0 )
      : Propagator(home), y(y0), x(x0), D(D0), n(p.n), m(p.m), k(p.k) {
         x.subscribe(home,*this,PC_INT_DOM);
      }

   forceinline
      MultiPack::MultiPack(Space& home, bool shared, Pack& p) 
      : Propagator(home,shared, p), n(p.n), m(p.m), k(p.k) {
         y.update(home, shared, p.y);
         x.update(home, shared, p.x);
         D.update(home, shared, p.D);
      }

   forceinline size_t 
      MultiPack::dispose(Space& home) {
         x.cancel(home,*this,PC_INT_DOM);
         (void) Propagator::dispose(home);
         return sizeof(*this);
      }
}}}

// STATISTICS: int-prop


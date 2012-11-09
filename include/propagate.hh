namespace Gecode { namespace Int { namespace CostBinPacking {

   /*
    * Item
    *
    */
   forceinline
      Item::Item(void) 
      : s(0) {}
   forceinline
      Item::Item(IntView b, int s0)
      : DerivedView<IntView>(b), s(s0) {}

   forceinline IntView 
      Item::bin(void) const {
         return x;
      }
   forceinline
      void Item::bin(IntView b) {
         x = b;
      }
   forceinline int 
      Item::size(void) const {
         return s;
      }
   forceinline void 
      Item::size(int s0) {
         s = s0;
      }

   forceinline void
      Item::update(Space& home, bool share, Item& i) {
         x.update(home,share,i.x);
         s = i.s;
      }


   forceinline bool 
      same(const Item& i, const Item& j) {
         return same(i.bin(),j.bin()) && (i.size() == j.size());
      }
   forceinline bool
      before(const Item& i, const Item& j) {
         return before(i.bin(),j.bin())
            || (same(i.bin(),j.bin()) && (i.size() == j.size()));
      }

   /// For sorting according to size
   forceinline bool 
      operator <(const Item& i, const Item& j) {
         return i.size() > j.size();
      }


   /*
    * Size set
    *
    */
   forceinline
      SizeSet::SizeSet(void) {}
   forceinline
      SizeSet::SizeSet(Region& region, int n_max) 
      : n(0), t(0), s(region.alloc<int>(n_max)) {}
   forceinline void
      SizeSet::add(int s0) {
         t += s0; s[n++] = s0;
      }
   forceinline int
      SizeSet::card(void) const {
         return n;
      }
   forceinline int
      SizeSet::total(void) const {
         return t;
      }
   forceinline int
      SizeSet::operator [](int i) const {
         return s[i];
      }

   forceinline
      SizeSetMinusOne::SizeSetMinusOne(void) {}
   forceinline
      SizeSetMinusOne::SizeSetMinusOne(Region& region, int n_max) 
      : SizeSet(region,n_max), p(-1) {}
   forceinline void
      SizeSetMinusOne::minus(int s0) {
         // This rests on the fact that items are removed in order
         do
            p++;
         while (s[p] > s0);
         assert(p < n);
      }
   forceinline int
      SizeSetMinusOne::card(void) const {
         assert(p >= 0);
         return n - 1;
      }
   forceinline int
      SizeSetMinusOne::total(void) const {
         assert(p >= 0);
         return t - s[p];
      }
   forceinline int
      SizeSetMinusOne::operator [](int i) const {
         assert(p >= 0);
         return s[(i < p) ? i : i+1];
      }



   /*
    * Packing propagator
    *
    */

   forceinline
      Pack::Pack(Home home, ViewArray<OffsetView>& l0, ViewArray<Item>& x0, IntView& z0)
      : Propagator(home), l(l0), x(x0), z(z0), t(0) {
         //l.subscribe(home,*this,PC_INT_BND);
         x.subscribe(home,*this,PC_INT_BND);
         for (int i=x.size(); i--; )
            t += x[i].size();
      }

   forceinline
      Pack::Pack(Space& home, bool shared, Pack& p) 
      : Propagator(home,shared,p), t(p.t) {
         l.update(home,shared,p.l);
         x.update(home,shared,p.x);
         z.update(home,shared,p.z);
      }

   forceinline size_t 
      Pack::dispose(Space& home) {
         //l.cancel(home,*this,PC_INT_BND);
         x.cancel(home,*this,PC_INT_BND);
         (void) Propagator::dispose(home);
         return sizeof(*this);
      }

}}}

// STATISTICS: int-prop


#include <cost_binpacking.hh>

namespace Gecode {

  void
  cost_binpacking(Home home, 
             const IntVarArgs& l, 
             const IntVarArgs& b, 
             const IntVar& z,
             const IntArgs& s,
             const IntSharedArray& D,
             IntConLevel) {
    using namespace Int;
    if (l.same(home,b))
      throw ArgumentSame("Int::cost_binpacking");
    if (b.size() != s.size())
      throw ArgumentSizeMismatch("Int::cost_binpacking");      
    if (b.size() * l.size() != D.size() )
      throw ArgumentSizeMismatch("Int::cost_binpacking");      

    for (int i=s.size(); i--; )
      Limits::positive(s[i],"Int::binpacking");
    
    if (home.failed()) return;

    ViewArray<OffsetView> lv(home,l.size());
    for (int i=l.size(); i--; )
      lv[i] = OffsetView(l[i],0);

    ViewArray<CostBinPacking::Item> bs(home,b.size());
    for (int i=bs.size(); i--; )
      bs[i] = CostBinPacking::Item(b[i],s[i]);

    //Support::quicksort(&bs[0], bs.size());

    IntView zv(z);
    GECODE_ES_FAIL(Int::CostBinPacking::Pack::post(home,lv,bs,zv,D));
  }
}

// STATISTICS: int-post

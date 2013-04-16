#include <cost_binpacking.hh>

namespace Gecode {
   void cost_binpacking(Home home, 
         int n, int m, int k,
         const IntVarArgs& y, 
         const IntVarArgs& x, 
         const IntSharedArray& D,
         IntConLevel) 
   {
      using namespace Int;

      if (n*k != D.size() )
         throw ArgumentSizeMismatch("Int::cost_multibin");      

      if (home.failed()) return;

      ViewArray<IntView> yv(home,y.size());
      for (int i=y.size(); i--; )
         yv[i] = IntView(y[i]);

      ViewArray<IntView> xv(home,x.size());
      for (int i=xv.size(); i--; )
         xv[i] = IntView(x[i]);

      GECODE_ES_FAIL(Int::CostBinPacking::MultiPack::post(home,n,m,k,yv,xv,D));
   }
}

// STATISTICS: int-post

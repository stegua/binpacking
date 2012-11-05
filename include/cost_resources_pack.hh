#ifndef __MY_COST_RESOURCES_
#define __MY_COST_RESOURCES_

#include "my_types.hh"

///--------------------------------------------------------------------------------
/// Resource and cost data
struct CostResources {
   cost_t     c;
   resources  r;

   CostResources(void) {}

   void setData( cost_t _c, int k ) {
      c = _c; 
      r.resize(k,0);
   }

   void add( const CostResources& o, const Arc& a ) {
      assert( r.size() == o.r.size() );

      c = o.c + a.c;
      int k = r.size();
      for ( int i = 0; i < k; ++i ) 
         r[i] = o.r[i] + a.r[i];
   }

   inline cost_t estimateCost ( const Arc& a, const CostResources& b ) const {
      return (c + a.c + b.c);
   }

   inline cost_t estimateCost ( const Arc& a, cost_t b ) const {
      return (c + a.c + b);
   }

   inline cost_t computeCost ( const Arc& a ) const {
      return (c + a.c);
   }

   bool isPathFeasible ( const Arc& a, const resources& B, const resources& U) const {
      int k = r.size();
      for ( int l = 0; l < k; ++l )
         if ( r[l] + a.r[l] + B[l] > U[l] )
            return false;
      return true;
   }

   bool isPathFeasible ( const Arc& a, const Arc& b, const CostResources& B, const resources& U) const {
      int k = r.size();
      for ( int l = 0; l < k; ++l )
         if ( r[l] + a.r[l] + b.r[l] + B.r[l] > U[l] )
            return false;
      return true;
   }
   bool isPathFeasible ( const Arc& a, const CostResources& B, const resources& U) const {
      int k = r.size();
      for ( int l = 0; l < k; ++l )
         if ( r[l] + a.r[l] + B.r[l] > U[l] )
            return false;
      return true;
   }
};


inline cost_t computeCost( cost_t c ) {
   return (c);
}

#endif /// __MY_COST_RESOURCES_

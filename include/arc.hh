#ifndef __MY_ARC_DIGRAPH_
#define __MY_ARC_DIGRAPH_

#include "my_types.hh"

using std::pair;

#include <boost/intrusive/list.hpp>
typedef boost::intrusive::list_member_hook< boost::intrusive::link_mode<boost::intrusive::normal_link> >  MyHook;

///------------------------------------------------------------------------------------------
/// My arc class
class Arc {
   public:
      /// Hooks ofr intrusive containers
      MyHook fs_hook;
      MyHook bs_hook;

      edge_t    id; /// Unique ID of the arc (to use as an hash)
      node_t    v;  /// Source node
      node_t    w;  /// Target node
      cost_t    c;  /// Cost of the arc
      cost_t    d;  /// Reduced cost of the arc (optimal multiplier)
      resources r;  /// Vector of resources consumed on the arc

      Arc(void) {}

      void setData ( edge_t _id, node_t _v, node_t _w, cost_t _c, const resources& _r ) {
         id = _id;  v = _v;  w = _w;  c = _c;  d = 0.0; r = _r;
      }

      /// TO be used to order the arc by their cost
      friend bool operator <(const Arc &l, const Arc &r)  {  return l.d < r.d;  }
      friend bool operator >(const Arc &l, const Arc &r)  {  return l.d > r.d;  }
};

/// Forward and Backward star: intrusive list
typedef boost::intrusive::list<Arc, boost::intrusive::member_hook<Arc, MyHook, &Arc::fs_hook > > FSArcList;
typedef boost::intrusive::list<Arc, boost::intrusive::member_hook<Arc, MyHook, &Arc::bs_hook > > BSArcList;
typedef FSArcList::iterator              FSArcIter;
typedef BSArcList::iterator              BSArcIter;
typedef std::pair<FSArcIter, FSArcIter>  FSArcIterPair;
typedef std::pair<BSArcIter, BSArcIter>  BSArcIterPair;


///--------------------------------------------------------------------------------
/// Methods for accesing different elements of an arc
#include <functional>
using std::unary_function;

template<class T>
struct CostView : public unary_function<const T&, cost_t> {
   inline cost_t operator()(const T& a) const { return a.c; }
};

template<class T>
struct GradientView : public unary_function<const T&, cost_t> {
   inline cost_t operator()(const T& a) const { return a.d; }
};

template<class T>
struct ResView : public unary_function<const T&, resource_t> {
   public:
      explicit ResView(int _k) : k(_k) {}
      inline resource_t operator()(const T& a) const { return a.r[k]; }  
   private:
      const int k;
};

template<class T>
struct NegResView : public unary_function<const T&, resource_t> {
   public:
      explicit NegResView(int _k) : k(_k) {}
      inline resource_t operator()(const T& a) const { return -a.r[k]; }  
   private:
      const int k;
};


/// Other typedefs
typedef CostView<Arc>              ArcCostView;
typedef GradientView<Arc>          ArcGradView;
typedef ResView<Arc>               ArcResView;
typedef NegResView<Arc>            NegArcResView;

#endif /// __MY_ARC_DIGRAPH_

#ifndef __MY_NODE_
#define __MY_NODE_

#include "arc.hh"

using std::pair;
using std::make_pair;

#include <boost/intrusive/list.hpp>
typedef boost::intrusive::list_member_hook< boost::intrusive::link_mode<boost::intrusive::normal_link> >  MyHook;

///------------------------------------------------------------------------------------------
/// My node class
class Node {
   public:
      /// Intrusive list hook
      MyHook n_hook;
      /// Unique ID of the node_t (to use as an hash)
      node_t id;       

      /// Standard constructors
      Node( void ) {}
      /// Copy constructor: OCCHIO NON COPI FS e BS
      Node( const Node& o ) : n_hook(o.n_hook), id(o.id) {}
      /// Assignment operator
      Node& operator=( const Node& other )
      {
         if( this == &other )
            return *this;
         this->~Node();
         new( this ) Node( other );
         return *this;
      } 
      /// Set the data
      inline void setData( node_t _id ) { id = _id; }
      /// Equality operator
      inline bool operator==(const Node& o) const { return (id == o.id);}  /// CONFRONTO ANCHE LE STAR?
      /// Basic getter
      inline node_t out_degree(void) const { return FS.size(); }
      inline node_t  in_degree(void) const { return BS.size(); }
      /// Iterators
      inline FSArcIterPair getIterFS(void) { return make_pair(FS.begin(), FS.end()); }  
      inline BSArcIterPair getIterBS(void) { return make_pair(BS.begin(), BS.end()); }  
      // /// Setters
      void addForwArc( Arc& a ) { FS.push_back(a); }
      void addBackArc( Arc& a ) { BS.push_back(a); }

      void removeForwArc( FSArcIter a ) { FS.erase( a ); }

      void removeForwArc( Arc& a ) { FS.erase( FS.iterator_to(a) ); }
      void removeBackArc( Arc& a ) { BS.erase( BS.iterator_to(a) ); }

      /// Sort by (REDUCED!) costs
      void sort( void ) { FS.sort(); }

   private:
      FSArcList FS;  /// Forward star  (set of arcs)
      BSArcList BS;  /// Backward star (set of arcs)
};

/// Node intrusive list
typedef boost::intrusive::list<Node, boost::intrusive::member_hook<Node, MyHook, &Node::n_hook > > NodeList;
typedef NodeList::iterator              NodeIter;
typedef std::pair<NodeIter, NodeIter>   NodeIterPair;

#endif /// __MY_NODE_

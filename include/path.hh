#ifndef __RCSP_PATH_
#define __RCSP_PATH_

#include "arc.hh"
#include "my_types.hh"
#include "dag_pack.hh"

#include <list>
using std::list;

#include <vector>
using std::vector;

class Path {
public:
   /// Cost of the path
   cost_t c; 
   /// Reduced Cost of the path
   cost_t d; 
   /// Resource consumption
   resources Rc;
   /// List of nodes in the path
   list< node_t >  path;
   /// Default constructor during filtering
//   Path ( node_t Source, node_t Target, cost_t _c, const Arc& a, const vector< node_t >& Pf, const vector< node_t >& Pb );
   /// Path costructore standard
   Path ( const DAG& G, node_t Source, node_t Target, const vector< edge_t >& Pf );
   /// Path constructor for enumeration
  // Path ( node_t Target, cost_t _c, const stack<node_t>& theStack);
   /// Empty constructor
   explicit Path ( void ) : c(0), d(0) {}
   /// Costruttore di copia
   Path ( const Path& p ) : c(p.c), d(p.d), Rc(p.Rc), path(p.path) {}
   /// Operatore di assegnamento
   Path& operator=( const Path& other ) {
      if( this == &other )
         return *this;
      this->~Path();
      new( this ) Path( other );
      return *this;
   }
};

//typedef ptr_vector<Path>  Paths;

#endif  /// close the __RCSP_PATH_


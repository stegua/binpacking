#ifndef __RCSP_PATH_
#define __RCSP_PATH_

#include "arc.hh"
#include "my_types.hh"
//#include "dag_pack.hh"
class DAG;
#include <vector>
using std::vector;

#include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

class Path {
public:
   /// Cost of the path
   cost_t c; 
   /// Reduced Cost of the path
   cost_t d; 
   /// Resource consumption
   resources Rc;
   /// List of nodes in the path
   vector< edge_t >  path;
   /// Path costructore standard
   Path ( const DAG& G, node_t Source, node_t Target, const vector< edge_t >& Pf );
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
   /// Update path cost with new weights
   cost_t updateCost( const DAG& G );
};

typedef ptr_vector<Path>  Paths;

#endif  /// close the __RCSP_PATH_


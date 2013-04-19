#ifndef __MY_TYPES_
#define __MY_TYPES_

/// My typedefs
#include <boost/cstdint.hpp>
#include <boost/integer_traits.hpp>

typedef int16_t       node_t;
typedef int16_t      edge_t;
typedef double    cost_t;
typedef int16_t      resource_t;
typedef double    dist_t;

#include <vector>
typedef std::vector<resource_t> resources;

#endif /// __MY_TYEPS_

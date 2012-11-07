#ifndef __MY_TYPES_
#define __MY_TYPES_

/// My typedefs
#include <boost/cstdint.hpp>
#include <boost/integer_traits.hpp>

typedef int32_t   node_t;
typedef int64_t   edge_t;
typedef float     cost_t;
typedef int32_t   resource_t;
typedef float     dist_t;

#include <vector>
typedef std::vector<resource_t> resources;

#endif /// __MY_TYEPS_

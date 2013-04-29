/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Christian Schulte, 2010
 *
 *  Last modified:
 *     $Date: 2010-10-06 23:20:35 +0200 (Wed, 06 Oct 2010) $ by $Author: schulte $
 *     $Revision: 11468 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef __GECODE_INT_MULTIBIN_PACKING_HH__
#define __GECODE_INT_MULTIBIN_PACKING_HH__

#include <gecode/int.hh>

namespace Gecode { namespace Int {
/**
 * \defgroup TaskModelIntBinPacking Multi Bin packing constraints
 * \ingroup TaskModelInt
 *
 * Constraint decomposition for modeling multidimensional bin packing problems. Propagation follows:
 *   Stefano Gualandi and Michele Lombardi. A Simple and Effective Decomposition for
 *   Multidimensional Binpacking Constraints. Tech. Report 2013.
 */
/* \brief Post propagator for multibin packing
 *
 * The constraint deals with \a n items, \a m bins, and \a k dimensions.
 * The variables in \a y are the loads for each bin and for each dimension, whereas the
 * variables in \a x define for each item into which bin it is packed.
 * The integer values \a s define the size of the items.
 *
 * It is propagated that for each \f$j\f$ with \f$0\leq j<|l|\f$ the
 * constraint \f$l_j=\sum_{0\leq i<|b|\wedge b_i=j}s_i\f$ holds and that
 * for each \f$i\f$ with \f$0\leq i<|b|\f$ the constraint
 * \f$0\leq b_i<|l|\f$ holds.
 *
 * Throws the following exceptions:
 *  - Of type Int::ArgumentSizeMismatch if \a b and \a s are not of
 *    the same size.
 *  - Of type Int::ArgumentSame if \a l and \a b share unassigned variables.
 *  - Of type Int::OutOfLimits if \a s contains a negative number.
 * 
 * \ingroup TaskModelIntBinPacking
 */
GECODE_INT_EXPORT void
multibinpacking(Home home,
      int n, int m, int k, 
      const IntVarArgs&      y, 
      const IntVarArgs&      x, 
      const IntSharedArray&  D,
      const IntSharedArray&  B,
      IntConLevel icl=ICL_DEF);
}}

#endif

// STATISTICS: int-prop


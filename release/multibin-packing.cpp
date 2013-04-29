/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Stefano Gualandi <stefano.gualandi@gmail.com>
 *
 *  Copyright:
 *     Christian Schulte, 2010
 *
 *  Last modified:
 *     $Date: 2012-09-07 17:31:22 +0200 (Fri, 07 Sep 2012) $ by $Author: schulte $
 *     $Revision: 13068 $
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

#include <gecode/int/distinct.hh>
#include <gecode/int/bin-packing.hh>

#include <multibin-packing.hh>

/// Use Cliquer to find all maximal clique in the conflict graph
extern "C" {
#include "cliquer.h"
#include "graph.h"
   /* Records a clique into the clique list using dynamic allocation.
    * Used as opts->user_function.
    */
   static int clique_count=0;
   static set_t *clique_list;
   static int clique_list_size=0;
   boolean record_clique_func(set_t s,graph_t *g,clique_options *opts) {
      if (clique_count>=clique_list_size) {
         clique_list=(set_t*)realloc(clique_list,(clique_list_size+512) * 
               sizeof(set_t));
         clique_list_size+=512;
      }
      clique_list[clique_count]=set_duplicate(s);
      clique_count++;
      return TRUE;
   }
}


namespace Gecode { namespace Int {

  void
  multibinpacking(Home home, 
         int n, int m, int k,
         const IntVarArgs& y, 
         const IntVarArgs& x, 
         const IntSharedArray& D,
         const IntSharedArray& B,
         IntConLevel) 
   {
      /// Check input sizes
      if (n*k != D.size() )
         throw ArgumentSizeMismatch("Int::multibinpacking");      

      if (k != B.size() )
         throw ArgumentSizeMismatch("Int::multibinpacking");      
    
      if (n != x.size() )
         throw ArgumentSizeMismatch("Int::multibinpacking");      
      
      if (m*k != y.size() )
         throw ArgumentSizeMismatch("Int::multibinpacking");      
      
      for (int i=B.size(); i--; )
         Limits::nonnegative(B[i],"Int::multibinpacking");

      if (home.failed()) return;

      /// Post first each single binpacking constraint
      /// Capacity constraint for each dimension
      for ( int j = 0; j < m; ++j )
         for ( int l = 0; l < k; ++l ) {
            IntView yi(y[j*k+l]);
            GECODE_ME_FAIL(yi.lq(home,B[l]));
         }

      /// Post a binpacking constraints for each dimension
      for ( int l = 0; l < k; ++l ) {
         ViewArray<OffsetView> yv(home,m);
         for (int j=m; j--; )
            yv[j] = OffsetView(y[j*k+l],0);

         ViewArray<BinPacking::Item> xs(home,x.size());
         for (int i=xs.size(); i--; )
            xs[i] = BinPacking::Item(x[i],D[i*k+l]);

         Support::quicksort(&xs[0], xs.size());

         GECODE_ES_FAIL(Int::BinPacking::Pack::post(home,yv,xs));
      }

      /// Clique Finding and Alldifferent posting
      {
         /// Firt construct the conflict graph
         graph_t* g = graph_new(n);
         for ( int a = 0; a < n-1; ++a ) {
            for ( int b = a+1; b < n; ++b ) {
               int v = a;  /// The variable with smaller domain
               int w = b;
               unsigned int nl = 0;
               if ( x[a].size() > x[b].size() ) {
                  v = b;
                  w = a;
               }
               IntVarValues i(x[v]);
               IntVarValues j(x[w]);
               while ( i() ) {
                  if ( i.val() != j.val() ) {
                     if ( i.val() < j.val() )
                        break;
                     else
                        ++i;
                  } else {
                     for ( int l = 0; l < k; ++l )
                        if ( D[a*k+l] + D[b*k+l] > B[l] ) {
                           nl++;
                           break;
                        }
                     ++i; ++j;
                  }
               }
               if ( nl >= x[v].size() )
                  GRAPH_ADD_EDGE(g,a,b);
            }
         }
         /// Consitency cheking: look for the maximum clique
         clique_options* opts;
         opts = (clique_options*) malloc (sizeof(clique_options));
         opts->time_function=NULL;
         opts->reorder_function=reorder_by_default;
         opts->reorder_map=NULL;
         opts->user_function=NULL;
         opts->user_data=NULL;
         opts->clique_list=NULL;
         opts->clique_list_length=0;
         set_t s;
         s = clique_find_single ( g, 0, 0, TRUE, opts);
         if ( s != NULL ) {
            if ( set_size(s) > m ) {
                set_free(s);
                free(opts);
                graph_free(g);
                GECODE_ES_FAIL(ES_FAILED);
            }
            if ( true ) { //option == 1 ) {  FIXING
               for ( int a = 0, j = 0; a < n; ++a ) {
                  if ( SET_CONTAINS_FAST(s,a) ) {
                     assert( x[a].in(j) );
                     IntView xi(x[a]);
                     GECODE_ME_FAIL(xi.eq(home,j++));
                  }
               }
            }
         }
         if ( s!=NULL )
            set_free(s);
         /// List every maximal clique in the conflict graph
         opts->user_function=record_clique_func;
         clique_find_all ( g, 2, 0, TRUE, opts);
         for ( int c = 0; c < clique_count; c++ ) {
            ViewArray<IntView> xv(home, set_size(clique_list[c]));
            for ( int a = 0, idx = 0; a < n; ++a )
               if ( SET_CONTAINS_FAST(clique_list[c],a) )
                  xv[idx++] = x[a];
            GECODE_ES_FAIL(Distinct::Dom<IntView>::post(home,xv));
            set_free(clique_list[c]);
         }
         free(opts);
         graph_free(g);
      }
   }
}}

// STATISTICS: int-post

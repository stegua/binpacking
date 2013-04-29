#include "cost_multibin.hh"
#include <stdio.h>
#include <dag_pack.hh>
#include <boost/unordered_map.hpp>
using boost::unordered_map;

#include <boost/timer.hpp>

#include <boost/unordered_set.hpp>
using boost::unordered_set;

extern "C" {
#include "cliquer.h"
#include "graph.h"
 /* Records a clique into the clique list using dynamic allocation.
 *  * Used as opts->user_function.
 *   */
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


#include <gecode/int/distinct.hh>
namespace Gecode { namespace Int { namespace CostMultiBinPacking {

   PropCost MultiPack::cost(const Space&, const ModEventDelta&) const {
      return PropCost::crazy(PropCost::HI,x.size());
   }

   Actor* MultiPack::copy(Space& home, bool share) {
      return new (home) MultiPack(home,share,*this);
   }

  ExecStatus 
  MultiPack::propagate(Space& home, const ModEventDelta& med) {
    /// Check if SUBSUMED
    {
       int i = 0;
       while ( i<n && x[i].assigned()) { i++; };
       if ( i == x.size() )
          return home.ES_SUBSUMED(*this);
    }
    //return ES_FIX;
    /// check the ammisibility
    int K = k*m;
    vector<bool> BO(m, false);
    vector<int>  CR(m*k,0);
    int fixed = 0;
    for ( int i = 0; i < n; ++i ) {
       if ( x[i].assigned() ) {
          fixed++;
          int j = x[i].val();
          BO[j] = true;
          for ( int l = 0; l < m; ++l )
             CR[j*k+l] += D[i*k+l];
       }
    }

    //if ( fixed == 0 )
      // return ES_NOFIX;

    /// Clique detection
       graph_t* g = graph_new(n);
       int max_dom = 0;
       for ( int a = 0; a < n-1; ++a ) {
          if ( !x[a].assigned() ) {
             //g->weights[a] = 100+m-x[a].size();
             if ( x[a].size() > max_dom )
                max_dom = x[a].size();
             for ( int b = 0; b < n; ++b ) {
                if ( a != b && !x[b].assigned() && x[a].size() <= x[b].size() ) {
                   int nl = 0;
                   IntVarValues i(x[a]);
                   IntVarValues j(x[b]);
                   while ( i() ) {
                      if ( i.val() != j.val() ) {
                         if ( i.val() < j.val() )
                            break;
                         else
                            ++i;
                      } else {
                         for ( int l = 0; l < k; ++l )
                            if ( D[a*k+l] + D[b*k+l] > B[l] - CR[i.val()*k+l] ) {
                               nl++;
                               break;
                            }
                         ++i; ++j;
                      }
                   }
                   if ( nl >= x[a].size() )
                      GRAPH_ADD_EDGE(g,a,b);
                }
             }
          }
       }

       clique_options* opts;
       set_t s;
       opts = (clique_options*) malloc (sizeof(clique_options));
       opts->time_function=NULL;
       opts->reorder_function=reorder_by_default;
       opts->reorder_map=NULL;
       opts->user_function=NULL;
       opts->user_data=NULL;
       opts->clique_list=NULL;
       opts->clique_list_length=0;
       for ( int a = 0; a < n; ++ a ) {
          if ( !x[a].assigned() ) {
             g->weights[a] = 100;
             s = clique_find_single ( g, 0, 0, TRUE, opts);
             if ( s != NULL && set_size(s) > 0 ) {
                max_dom = 0;
          //      fprintf(stdout,"clique %d\t", set_size(s));
                ViewArray<IntView> xv(home, set_size(s));
                int idx = 0;
                for ( int a = 0; a < n; ++a )
                   if ( SET_CONTAINS_FAST(s,a) ) {
                      xv[idx] = x[a];
            //          printf("x[%d]=%d ", a, x[a].size());
                      if ( x[a].size() > max_dom )
                         max_dom = x[a].size();
                      idx++;
                   }
                /*ExecStatus es = Distinct::prop_bnd<IntView>(home,xv);
                GECODE_ES_CHECK(es);
                if ( xv.size() > 2 )
                   es = Distinct::prop_val<IntView,true>(home,xv);
                GECODE_ES_CHECK(es);
                printf("\n");*/
             }
             if ( s != NULL && set_size(s) > max_dom ) {
                set_free(s);
                free(opts);
                graph_free(g);
                fprintf(stdout,"Failed\n");
                return ES_FAILED;
             }
             g->weights[a] = 1;
          }
       }
       if ( s!=NULL )
          set_free(s);
       free(opts);
       graph_free(g);

       return ES_NOFIX;

       /// Create the graph and propagates: x = x
       int N = n*m+2;
       int M = (n-1)*(m*m)+2*m;
       int S = N-2;
       int T = N-1;

    resources U(K,0);
    for ( int l = 0; l < K; ++l )
       U[l] = y[l].max();

    DAG G (N, M, U);
    
    /// Build the Directed Acyclic Graph
    typedef std::pair<int, int> NodePair;
    unordered_map< NodePair, edge_t > As;
    for ( int i = 0; i < n-1; ++i ) {
       for ( IntVarValues j(x[i+1]); j(); ++j ) {
          resources R(K,0);
          for ( int l = 0; l < k; ++l )
             R[j.val()*k+l] = D[(i+1)*k+l];//A[i+1][l];
          for ( IntVarValues h(x[i]); h(); ++h ) {
             edge_t arc_id = G.addArc( i*m+h.val(), (i+1)*m+j.val(), 0, R );
             As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))] = arc_id;
          }
       }
    }
    /// Arcs from the source node
    for ( IntVarValues j(x[0]); j(); ++j ) {
       resources R(K,0);
       for ( int l = 0; l < k; ++l )
          R[j.val()*k+l] = D[l];//x[0].size();
       edge_t arc_id = G.addArc( S, j.val(), 0, R );
       As[make_pair(S, j.val())] = arc_id;
    }
    /// Arcs to the destination node
    for ( IntVarValues j(x[n-1]); j(); ++j ) {
       resources R(K,0);
       G.addArc( (n-1)*m+j.val(), T, 0, R );
    }

    boost::timer TTIMER;
    double filter_time = 0;
    /// Filter the arcs
    cost_t LB;
    cost_t UB;
    int status = 0;
    for ( int q = 0; q < m; ++q ) {
       bool flag = false;
       if ( BO[q] ) {
          if ( !G.isLPdefined() )
             G.setLPmodel(S, T);
          for ( int l = 0; l < k; ++l ) {
             if ( !y[q*k+l].assigned() && y[q*k+l].min() > 0 ) {
                LB = 0;
                UB = y[q*k+l].max();
                //fprintf(stdout,"%.1f ", UB);
                /// Archi da item a item
                for ( int i = 0; i < n-1; ++i ) {
                   for ( IntVarValues j(x[i+1]); j(); ++j ) {
                      for ( IntVarValues h(x[i]); h(); ++h ) {
                         if ( j.val() == q )
                            G.setArcCost(As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))], D[(i+1)*k+l]);
                         else
                            G.setArcCost(As[make_pair(i*m+h.val(), ((i+1)*m+j.val()))], 0);
                      }
                   }
                }
                /// Archi dalla sorgente 
                for ( IntVarValues j(x[0]); j(); ++j ) {
                   if ( j.val() == q )
                      G.setArcCost( As[make_pair(S, j.val())], D[l]);
                   else
                      G.setArcCost( As[make_pair(S, j.val())], 0);
                }

                double t0 = TTIMER.elapsed();
                status = G.filter(S,T,LB,UB);
                filter_time += TTIMER.elapsed()-t0;

                if ( status == 2 ) {
                   fprintf(stdout,"fail-reach\n");
                   return ES_FAILED;
                }

                /// Aumenta il lower bound della variable di load
                if ( status == 0 ) {
                   int LB0 = int(LB)-1;
                   if ( fabs(LB-ceil(LB)) > EPS )
                      LB0 = int(ceil(LB))-1;
                   else 
                      LB0 = int(round(LB))-1;
                   if ( y[q*k+l].max() < LB0 ) {
                      //fprintf(stdout,"fail-bound\n");
                      return ES_FAILED;
                   }
                   if ( y[q*k+l].min() < LB0 ) {
                      //fprintf(stdout,"prune #%.1f %d %d\n", LB, LB0, y[q*k+l].min());
                      GECODE_ME_CHECK(y[q*k+l].gq(home,LB0));
                      flag = true;
                   }
                }
             }
             //fprintf(stdout,"\n");
          }
       }
       //if ( flag )
       // break;
    }

    if ( ES_FAILED == G.filterArcs(n,m,k,x,home) )
       return ES_FAILED;
    //fprintf(stdout, "Time %.5f %.5f\n",  TTIMER.elapsed(), filter_time);    
    return ES_FIX;
  }

  ExecStatus
     MultiPack::post(Home home, int n, int m, int k,
           ViewArray<IntView>& y, ViewArray<IntView>& x, 
           const IntSharedArray& D, const IntSharedArray& B) {
        if (x.size() == 0) {
           // No items to be packed
           for (int i=y.size(); i--; )
              GECODE_ME_CHECK(y[i].eq(home,0));
           return ES_OK;
        } else if (y.size() == 0) {
           // No bins available
           return ES_FAILED;
        } else {
           // Constrain bins 
           for (int i=x.size(); i--; ) {
              GECODE_ME_CHECK(x[i].gq(home,0));
              GECODE_ME_CHECK(x[i].le(home,y.size()));
           }
           {
              vector< unordered_set<int> > CLIQUES;
              graph_t* g = graph_new(n);
              int max_dom = 0;
              for ( int a = 0; a < n-1; ++a ) {
                 if ( true ) 
                    g->weights[a] = n;
                 if ( !x[a].assigned() ) {
                    if ( x[a].size() > max_dom )
                       max_dom = x[a].size();
                    for ( int b = 0; b < n; ++b ) {
                       if ( a != b && !x[b].assigned() && x[a].size() <= x[b].size() ) {
                          int nl = 0;
                          IntVarValues i(x[a]);
                          IntVarValues j(x[b]);
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
                          if ( nl >= x[a].size() )
                             GRAPH_ADD_EDGE(g,a,b);
                       }
                    }
                 }
              }

              clique_options* opts;
              set_t s;
              opts = (clique_options*) malloc (sizeof(clique_options));
              opts->time_function=NULL;
              opts->reorder_function=reorder_by_random;
              opts->reorder_map=NULL;
              opts->user_function=record_clique_func;
              opts->user_data=NULL;
              opts->clique_list=NULL;
              opts->clique_list_length=0;
              for ( int a = 0; a < n; ++a ) {
                 if ( true ) {
                    clique_find_all ( g, 2, 0, TRUE, opts);
                    for ( int c = 0; c < clique_count; c++ ) {
                       ViewArray<IntView> xv(home, set_size(clique_list[c]));
                       printf("Clique size %d\n", set_size(clique_list[c]));
                       int idx = 0;
                       for ( int a = 0; a < n; ++a )
                          if ( SET_CONTAINS_FAST(clique_list[c],a) ) {
                             xv[idx] = x[a];
                             idx++;
                          }
                       if ( set_size(clique_list[c]) > m )
                          return ES_FAILED;
                       else
                          GECODE_ES_CHECK(Distinct::Dom<IntView>::post(home,xv));
                    }
                    printf("Alldiff %d\n", clique_count);
                    break;
                 } else { /// one clique per ndoe
                    g->weights[a] += 100;
                    s = clique_find_single ( g, 0, 0, TRUE, opts);
                    if ( s != NULL && set_size(s) > 0 ) {
                       fprintf(stdout,"%d: clique %d\t", a, set_size(s));
                       ViewArray<IntView> xv(home, set_size(s));
                       int idx = 0;
                       unordered_set<int> cli;
                       for ( int a = 0; a < n; ++a )
                          if ( SET_CONTAINS_FAST(s,a) ) {
                             cli.insert(a);
                             g->weights[a] -= 1;
                             xv[idx] = x[a];
                             printf("%d ", a);
                             idx++;
                          }
                       if ( set_size(s) > m )
                          return ES_FAILED;
                       if ( set_size(s) > 1 ) {
                          bool flag = false;
                          for ( int t = 0; t < CLIQUES.size(); ++t )
                             if (CLIQUES[t] == cli) {
                                flag = true;
                                break;
                             }
                          if ( !flag ) {
                             CLIQUES.push_back(cli);
                             GECODE_ES_CHECK(Distinct::Dom<IntView>::post(home,xv));
                          }
                       }
                       printf("\n");
                    }
                    g->weights[a] -= 100;
                 }
              }
              if ( s!=NULL )
                 set_free(s);
              free(opts);
              graph_free(g);
              printf("Alldiff %d\n", CLIQUES.size());
           } 
           //      (void) new (home) MultiPack(home, n, m, k, y, x, D, B);
           return ES_OK;
        }
        /// Check also the z variable!
     }
}}}

// STATISTICS: int-prop


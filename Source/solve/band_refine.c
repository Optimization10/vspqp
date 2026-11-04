#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "core.h"
#include "solve_internal.h"
#include "utils.h"
#include "vs.h"


// Two-sided band refinement of a vertex separator (Sanders and Schulz).
//
// The band is the separator S plus one BFS into each shore. The X side stops
// for good at the first vertex that would make c(A cap X) exceed
// u - c(Y) - c(S), the Y side likewise at u - c(X) - c(S), and both stop once
// the total degree of the band would exceed the degree cap. The rest of X is
// contracted into the source and the rest of Y into the sink. Each band vertex
// v is split into v_in -> v_out with capacity c(v), and graph edges inside the
// band get infinite capacity. A minimum s-t cut, found by Edmonds-Karp or by
// igraph (band_flow_type), is then a minimum weight separator inside the
// band. Unlike a cover step, it can move into both shores at once, so it is
// not limited by the slack of one shore.
//
// The cut ignores the bound u on the shores. Of the two extreme minimum cuts
// (closest to the source and closest to the sink), the feasible one with the
// lightest separator, then the better balance, is taken. If neither is
// feasible, the pass makes no move. A move is accepted if the separator gets
// lighter, or stays the same weight and the shores get more balanced.

#define S_NODE  0
#define T_NODE  1
#define IN(li)  (2 + 2*(li))
#define OUT(li) (3 + 2*(li))

// 0 = X, 1 = Y, 2 = separator
#define SIDE(x, y, v) ((x)[v] ? 0 : ((y)[v] ? 1 : 2))

// The infinite capacity, 2 * total weight + 1, can exceed VSPLONG
typedef long long band_cap ;


/* ========================================================================== */
// Flow network: forward star, arcs stored in pairs e, e^1
/* ========================================================================== */

typedef struct
{
    int       nn, ne ;
    int      *head, *nxt, *to ;
    band_cap *cap ;   // residual capacities
    int      *it, *q ;
} Band_Net ;


static void Net_Init(Band_Net *N, int nn, long long max_ne)
{
    N->nn    = nn ;
    N->ne    = 0 ;
    N->head  = malloc(nn * sizeof(int)) ;
    N->it    = malloc(nn * sizeof(int)) ;
    N->q     = malloc(nn * sizeof(int)) ;
    N->nxt   = malloc(max_ne * sizeof(int)) ;
    N->to    = malloc(max_ne * sizeof(int)) ;
    N->cap   = malloc(max_ne * sizeof(band_cap)) ;
    for(int i = 0 ; i < nn ; i++) N->head[i] = -1 ;
}


static void Net_Free(Band_Net *N)
{
    free(N->head) ;
    free(N->it) ;
    free(N->q) ;
    free(N->nxt) ;
    free(N->to) ;
    free(N->cap) ;
}


// Arc u -> v with capacity c, and its reverse with capacity 0
static void Net_Add(Band_Net *N, int u, int v, band_cap c)
{
    int e = N->ne ;

    N->to[e]   = v ;
    N->cap[e]  = c ;
    N->nxt[e]  = N->head[u] ;
    N->head[u] = e ;

    N->to[e+1]   = u ;
    N->cap[e+1]  = 0 ;
    N->nxt[e+1]  = N->head[v] ;
    N->head[v]   = e + 1 ;

    N->ne += 2 ;
}


// Edmonds-Karp on the residual capacities in N->cap, as in
// Min_Weighted_Vtx_Cover (vs.c): one BFS per augmenting path, stopped as soon
// as t is discovered. N->it holds the arc by which each node was reached
static void Net_Max_Flow_EK(Band_Net *N, int s, int t)
{
    int *pe = N->it ;
    int *q  = N->q ;

    while(true)
    {
        /* === Find an augmenting path ====================================== */
        int head = 0, tail = 0 ;
        for(int i = 0 ; i < N->nn ; i++) pe[i] = -2 ;
        pe[s] = -1 ;
        q[tail++] = s ;
        while(head < tail && pe[t] == -2)
        {
            int u = q[head++] ;
            for(int e = N->head[u] ; e != -1 ; e = N->nxt[e])
            {
                int v = N->to[e] ;
                if(N->cap[e] == 0 || pe[v] != -2) continue ;
                pe[v] = e ;
                if(v == t) break ;
                q[tail++] = v ;
            }
        }
        if(pe[t] == -2) break ;

        /* === Augment by the bottleneck capacity =========================== */
        band_cap f = N->cap[pe[t]] ;
        for(int v = t ; v != s ; v = N->to[pe[v] ^ 1])
        {
            if(N->cap[pe[v]] < f) f = N->cap[pe[v]] ;
        }
        for(int v = t ; v != s ; v = N->to[pe[v] ^ 1])
        {
            N->cap[pe[v]]     -= f ;
            N->cap[pe[v] ^ 1] += f ;
        }
    }
}


// mark[v] = 1 if v is reachable from s in the residual graph
static void Net_Reach_From(Band_Net *N, int s, char *mark)
{
    int head = 0, tail = 0 ;
    memset(mark, 0, N->nn) ;
    mark[s] = 1 ;
    N->q[tail++] = s ;
    while(head < tail)
    {
        int u = N->q[head++] ;
        for(int e = N->head[u] ; e != -1 ; e = N->nxt[e])
        {
            int v = N->to[e] ;
            if(N->cap[e] > 0 && !mark[v])
            {
                mark[v] = 1 ;
                N->q[tail++] = v ;
            }
        }
    }
}


// mark[v] = 1 if t is reachable from v in the residual graph
static void Net_Reach_To(Band_Net *N, int t, char *mark)
{
    int head = 0, tail = 0 ;
    memset(mark, 0, N->nn) ;
    mark[t] = 1 ;
    N->q[tail++] = t ;
    while(head < tail)
    {
        int u = N->q[head++] ;
        for(int e = N->head[u] ; e != -1 ; e = N->nxt[e])
        {
            // e^1 is the arc to[e] -> u
            int v = N->to[e] ;
            if(N->cap[e ^ 1] > 0 && !mark[v])
            {
                mark[v] = 1 ;
                N->q[tail++] = v ;
            }
        }
    }
}


/* ========================================================================== */
// Logging: output/band_calls.csv (one row per Band_Refine call) and
// output/band_passes.csv (one row per band pass). Only writes files; the
// results are unchanged.
/* ========================================================================== */

static double Band_Now(void)
{
    struct timespec ts ;
    clock_gettime(CLOCK_MONOTONIC, &ts) ;
    return ts.tv_sec + ts.tv_nsec * 1e-9 ;
}

// Open path for appending; write header first if the file is new or empty
static FILE *Band_Log_Open(const char *path, const char *header)
{
    FILE *f = fopen(path, "a") ;
    if(f && ftell(f) == 0) fprintf(f, "%s\n", header) ;
    return f ;
}

static const char *Band_Problem_Name(void)
{
    return G_Current_Problem_Name ? G_Current_Problem_Name : "NA" ;
}

/* ========================================================================== */
// Band
/* ========================================================================== */

typedef struct
{
    Band_Net  N ;
    int       nb ;              // number of band vertices
    int      *band ;            // band[li] is a vertex of the graph
    int      *tmp ;             // sides of the band vertices for one cut
    int      *best ;            // sides of the band vertices for the best cut
    char     *zs, *zt ;
    VSPINT   *x, *y, *c ;
    VSPINT    u ;
    band_cap  wx0, wy0, sep0 ;  // current weights
    band_cap  best_sep, best_imb ;
    bool      have_best ;
    band_flow_type flow ;       // max flow algorithm
    double    t_flow ;          // logging: seconds in the last max flow
    band_cap  flow_val ;        // logging: value of the last max flow
} Band ;


// Sides of the band vertices for the cut whose source side is mark
// (invert = 0) or the complement of mark (invert = 1), and the new weights
static void Band_Eval
(
    Band       *B,
    const char *mark,
    int         invert,
    band_cap   *sep,
    band_cap   *wx,
    band_cap   *wy
)
{
    band_cap w[3] = {B->wx0, B->wy0, B->sep0} ;
    for(int li = 0 ; li < B->nb ; li++)
    {
        int v    = B->band[li] ;
        int zin  = mark[IN(li)]  ^ invert ;
        int zout = mark[OUT(li)] ^ invert ;
        int ns   = zin ? (zout ? 0 : 2) : 1 ;
        int os   = SIDE(B->x, B->y, v) ;
        if(ns != os)
        {
            w[os] -= B->c[v] ;
            w[ns] += B->c[v] ;
        }
        B->tmp[li] = ns ;
    }
    *wx  = w[0] ;
    *wy  = w[1] ;
    *sep = w[2] ;
}


// Max flow, then evaluate both extreme minimum cuts and keep the best feasible
// one if it improves on the best so far
static void Band_Solve(Band *B)
{
    Band_Net *N = &B->N ;

    double t0 = Band_Now() ;
    if(B->flow == BAND_IGRAPH)
    {
        Band_Igraph_Max_Flow(N->nn, N->ne, N->to, N->cap, S_NODE, T_NODE) ;
    }
    else Net_Max_Flow_EK(N, S_NODE, T_NODE) ;
    B->t_flow   = Band_Now() - t0 ;
    B->flow_val = 0 ;
    for(int e = N->head[S_NODE] ; e != -1 ; e = N->nxt[e]) B->flow_val += N->cap[e ^ 1] ;
    Net_Reach_From(N, S_NODE, B->zs) ;
    Net_Reach_To(N, T_NODE, B->zt) ;

    for(int ex = 0 ; ex < 2 ; ex++)
    {
        band_cap sep, wx, wy ;
        if(ex == 0) Band_Eval(B, B->zs, 0, &sep, &wx, &wy) ;  // nearest s
        else        Band_Eval(B, B->zt, 1, &sep, &wx, &wy) ;  // nearest t

        if(wx > B->u || wy > B->u) continue ;

        band_cap imb = (wx > wy) ? wx - wy : wy - wx ;
        if(sep < B->best_sep || (sep == B->best_sep && imb < B->best_imb))
        {
            B->best_sep  = sep ;
            B->best_imb  = imb ;
            B->have_best = true ;
            memcpy(B->best, B->tmp, B->nb * sizeof(int)) ;
        }
    }
}


// Move the band vertices to the sides in B->best. Az is updated column by
// column for the vertices that move
static void Band_Apply(Band *B, State *X, State *Y, Problem *P)
{
    VSPINT *c         = P->graph->c ;
    State  *states[2] = {X, Y} ;

    for(int li = 0 ; li < B->nb ; li++)
    {
        int v = B->band[li] ;
        for(int s = 0 ; s < 2 ; s++)
        {
            VSPINT *z    = states[s]->point->int_part ;
            VSPINT  want = (B->best[li] == s) ;
            if(z[v] == want) continue ;

            z[v] = want ;
            Az_pm_col_of_A(states[s]->Az->int_part, v, want ? PLUS : MINUS, P) ;
            states[s]->weight += want ? c[v] : -c[v] ;
        }
    }
}


#ifdef DEBUG
// Check that X, Y form a true separator: no vertex in both shores, no edge
// between the shores, and weights that agree. check_states covers Az and u.
static void Band_Check_Separator(State *X, State *Y, Problem *P)
{
    SET_FUNCTION_COUNTER ;

    VSPINT   i ;
    VSPLONG  p ;
    VSPINT   n      = P->graph->n ;
    VSPINT  *c      = P->graph->c ;
    VSPINT  *Ai     = P->graph->Ai ;
    VSPLONG *Ap     = P->graph->Ap ;
    VSPINT  *x      = X->point->int_part ;
    VSPINT  *y      = Y->point->int_part ;
    VSPINT   wsep   = 0 ;
    VSPINT   nboth  = 0 ;
    VSPINT   ncross = 0 ;

    check_states(X, Y, P) ;

    for(i = 0 ; i < n ; i++)
    {
        if(x[i] && y[i])   nboth++ ;
        if(!x[i] && !y[i]) wsep += c[i] ;
        if(!x[i]) continue ;

        for(p = Ap[i] ; p < Ap[i+1] ; p++)
        {
            if(y[Ai[p]]) ncross++ ;
        }
    }

    if(nboth > 0)  ERR("Vertex in both shores") ;
    if(ncross > 0) ERR("Edge between the shores") ;
    if(wsep != P->graph->total_vtx_cost - X->weight - Y->weight)
    {
        ERR("Separator weight disagrees with the shore weights") ;
    }
}
#endif


// Refine the separator X, Y with at most max_passes band moves. The weight
// bounds of each side limit the band, and its total degree (separator
// included) is at most deg_factor * sqrt(C(V)). If the separator alone
// exceeds that, or no move is accepted, Cover_Refine is run instead.
// Returns the number of band passes
int Band_Refine
(
    State   *X,
    State   *Y,
    Problem *P,
    VSPINT   max_bip,
    VSPINT   max_passes,
    band_flow_type flow,
    double   deg_factor
)
{
    SET_FUNCTION_COUNTER ;
    if(max_bip <= 0) return 0 ;

    // Logging
    double    lg_t0      = Band_Now() ;
    long long lg_sep_in  = (long long) P->graph->total_vtx_cost - X->weight - Y->weight ;
    long long lg_wx_in   = X->weight, lg_wy_in = Y->weight ;
    long long lg_degS    = -1 ;          // degree of the first separator
    int       lg_flows   = 0 ;           // passes that ran a max flow
    int       lg_acc     = 0 ;           // accepted passes
    double    lg_t_flows = 0 ;
    const char *lg_skip  = "" ;

    VSPINT   i ;
    VSPLONG  p ;
    VSPINT   n        = P->graph->n ;
    VSPINT  *c        = P->graph->c ;
    VSPINT  *Ai       = P->graph->Ai ;
    VSPLONG *Ap       = P->graph->Ap ;
    VSPINT  *x        = X->point->int_part ;
    VSPINT  *y        = Y->point->int_part ;
    VSPINT   tot      = P->graph->total_vtx_cost ;
    int      passes   = 0 ;
    bool     moved    = false ;
    // Cap on the total degree of the band, separator included
    long long deg_cap = (long long) (deg_factor * sqrt((double) tot)) ;

    VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;

    int nS = 0 ;
    for(i = 0 ; i < n ; i++) if(!x[i] && !y[i]) nS++ ;
    int lg_S_in = nS ;
    if(nS == 0)
    {
        double lg_tc = Band_Now() ;
        int lg_cov   = Cover_Refine(X, Y, P, max_bip) ;
        double lg_t1 = Band_Now() ;
        FILE *f = Band_Log_Open("output/band_calls.csv",
            "problem,level,n,u,wx_in,wy_in,sep_in,S_in,degS,deg_cap,skip,band_passes,band_accepted,"
            "sep_after_band,t_band,t_maxflow,cover_ran,cover_iters,cover_moves,sep_after_cover,t_cover,t_total") ;
        if(f)
        {
            long long sep_out = (long long) tot - X->weight - Y->weight ;
            fprintf(f, "%s,%lld,%lld,%lld,%lld,%lld,%lld,0,0,%lld,empty,0,0,%lld,0,0,1,%d,%d,%lld,%.6f,%.6f\n",
                Band_Problem_Name(), (long long) P->level, (long long) n, (long long) P->u,
                lg_wx_in, lg_wy_in, lg_sep_in, deg_cap, lg_sep_in,
                lg_cov, lg_cov > 0 ? lg_cov - 1 : 0, sep_out, lg_t1 - lg_tc, lg_t1 - lg_t0) ;
            fclose(f) ;
        }
        return 0 ;
    }

    Band B ;
    B.x     = x ;
    B.y     = y ;
    B.c     = c ;
    B.u     = P->u ;
    B.flow  = flow ;
    B.band  = malloc(n * sizeof(int)) ;
    B.tmp   = malloc(n * sizeof(int)) ;
    B.best  = malloc(n * sizeof(int)) ;
    B.zs    = malloc(2 * (size_t) n + 2) ;
    B.zt    = malloc(2 * (size_t) n + 2) ;
    int *Slist = malloc(n * sizeof(int)) ;
    int *loc   = malloc(n * sizeof(int)) ;   // band index, -1 if not in band
    for(i = 0 ; i < n ; i++) loc[i] = -1 ;

    nS = 0 ;
    for(i = 0 ; i < n ; i++) if(!x[i] && !y[i]) Slist[nS++] = i ;

    while(passes < max_passes)
    {
        passes++ ;

        /* === Band: S plus one BFS into each shore ======================= */
        int       nb      = 0, k ;
        long long deg_sum = 0 ;

        for(k = 0 ; k < nS ; k++)
        {
            loc[Slist[k]] = nb ;
            B.band[nb++]  = Slist[k] ;
            deg_sum      += Ap[Slist[k]+1] - Ap[Slist[k]] ;
        }
        if(passes == 1) lg_degS = deg_sum ;
        if(deg_sum > deg_cap)
        {
            if(passes == 1) lg_skip = "degS" ;
            for(k = 0 ; k < nb ; k++) loc[B.band[k]] = -1 ;
            break ;
        }

        // One BFS into X and one into Y, run together layer by layer (no X-Y
        // edges, so each side is its own BFS). The X side stops for good at
        // the first vertex that would make c(A cap X) exceed u - c(Y) - c(S),
        // the Y side likewise; a bound <= 0 gives that side nothing. Both stop
        // once the total degree would exceed the cap
        band_cap wS      = (band_cap) tot - X->weight - Y->weight ;
        band_cap bnd[2]  = {(band_cap) P->u - Y->weight - wS,
                            (band_cap) P->u - X->weight - wS} ;
        band_cap add[2]  = {0, 0} ;
        bool     stop[2] = {bnd[0] <= 0, bnd[1] <= 0} ;
        int lo = 0, hi = nb ;
        while(lo < hi && !(stop[0] && stop[1]))
        {
            for(k = lo ; k < hi && !(stop[0] && stop[1]) ; k++)
            {
                int v = B.band[k] ;
                for(p = Ap[v] ; p < Ap[v+1] ; p++)
                {
                    int w = Ai[p] ;
                    if(loc[w] != -1) continue ;

                    int sd = SIDE(x, y, w) ;
                    if(sd == 2 || stop[sd]) continue ;
                    if(add[sd] + c[w] > bnd[sd])
                    {
                        stop[sd] = true ;
                        continue ;
                    }

                    long long dw = Ap[w+1] - Ap[w] ;
                    if(deg_sum + dw > deg_cap)
                    {
                        stop[0] = stop[1] = true ;
                        break ;
                    }

                    loc[w]       = nb ;
                    B.band[nb++] = w ;
                    add[sd]     += c[w] ;
                    deg_sum     += dw ;
                }
            }
            lo = hi ;
            hi = nb ;
        }
        B.nb = nb ;

        /* === Network ======================================================= */
        long long max_ne = 3LL * nb ;
        for(k = 0 ; k < nb ; k++)
        {
            int v = B.band[k] ;
            for(p = Ap[v] ; p < Ap[v+1] ; p++)
            {
                if(Ai[p] != v && loc[Ai[p]] != -1) max_ne++ ;
            }
        }
        Net_Init(&B.N, 2 * nb + 2, 2 * max_ne) ;
        band_cap inf = 2LL * (band_cap) tot + 1 ;

        for(int li = 0 ; li < nb ; li++)
        {
            int v     = B.band[li] ;
            int to_s  = 0 ;
            int to_t  = 0 ;
            Net_Add(&B.N, IN(li), OUT(li), (band_cap) c[v]) ;
            for(p = Ap[v] ; p < Ap[v+1] ; p++)
            {
                int w = Ai[p] ;
                if(w == v) continue ;
                if(loc[w] != -1)
                {
                    Net_Add(&B.N, OUT(li), IN(loc[w]), inf) ;
                }
                else
                {
                    // Outside the band: X is the source and Y the sink. Every
                    // separator vertex is in the band
                    int sd = SIDE(x, y, w) ;
                    if(sd == 0) to_s = 1 ;
                    if(sd == 1) to_t = 1 ;
                }
            }
            if(to_s) Net_Add(&B.N, S_NODE, IN(li), inf) ;
            if(to_t) Net_Add(&B.N, OUT(li), T_NODE, inf) ;
        }

        /* === Solve ========================================================= */
        B.wx0       = X->weight ;
        B.wy0       = Y->weight ;
        B.sep0      = tot - X->weight - Y->weight ;
        B.best_sep  = B.sep0 ;
        B.best_imb  = (B.wx0 > B.wy0) ? B.wx0 - B.wy0 : B.wy0 - B.wx0 ;
        B.have_best = false ;

        Band_Solve(&B) ;
        lg_flows++ ;
        lg_t_flows += B.t_flow ;
        {
            FILE *f = Band_Log_Open("output/band_passes.csv",
                "problem,level,pass,S,A,degA,deg_cap,nodes,arcs,flow,t_maxflow,accepted,sep_before,sep_after") ;
            if(f)
            {
                fprintf(f, "%s,%lld,%d,%d,%d,%lld,%lld,%d,%d,%lld,%.6f,%d,%lld,%lld\n",
                    Band_Problem_Name(), (long long) P->level, passes, nS, nb, deg_sum, deg_cap,
                    B.N.nn, B.N.ne / 2, (long long) B.flow_val, B.t_flow, B.have_best ? 1 : 0,
                    (long long) B.sep0, (long long) (B.have_best ? B.best_sep : B.sep0)) ;
                fclose(f) ;
            }
        }

        Net_Free(&B.N) ;
        for(k = 0 ; k < nb ; k++) loc[B.band[k]] = -1 ;

        if(!B.have_best) break ;

        /* === Accept ======================================================== */
        Band_Apply(&B, X, Y, P) ;
        moved = true ;
        lg_acc++ ;
#ifdef DEBUG
        Band_Check_Separator(X, Y, P) ;
#endif

        // The new separator lies in the band
        nS = 0 ;
        for(k = 0 ; k < nb ; k++) if(B.best[k] == 2) Slist[nS++] = B.band[k] ;
        if(nS == 0) break ;
    }

    free(B.band) ;
    free(B.tmp) ;
    free(B.best) ;
    free(B.zs) ;
    free(B.zt) ;
    free(Slist) ;
    free(loc) ;

    double    lg_tb        = Band_Now() ;
    long long lg_sep_band  = (long long) tot - X->weight - Y->weight ;
    int       lg_cov       = -1 ;
    if(!moved) lg_cov = Cover_Refine(X, Y, P, max_bip) ;
    double    lg_t1        = Band_Now() ;
    {
        FILE *f = Band_Log_Open("output/band_calls.csv",
            "problem,level,n,u,wx_in,wy_in,sep_in,S_in,degS,deg_cap,skip,band_passes,band_accepted,"
            "sep_after_band,t_band,t_maxflow,cover_ran,cover_iters,cover_moves,sep_after_cover,t_cover,t_total") ;
        if(f)
        {
            long long sep_out = (long long) tot - X->weight - Y->weight ;
            fprintf(f, "%s,%lld,%lld,%lld,%lld,%lld,%lld,%d,%lld,%lld,%s,%d,%d,%lld,%.6f,%.6f,%d,%d,%d,%lld,%.6f,%.6f\n",
                Band_Problem_Name(), (long long) P->level, (long long) n, (long long) P->u,
                lg_wx_in, lg_wy_in, lg_sep_in, lg_S_in, lg_degS, deg_cap, lg_skip,
                lg_flows, lg_acc, lg_sep_band, lg_tb - lg_t0, lg_t_flows,
                lg_cov >= 0 ? 1 : 0, lg_cov >= 0 ? lg_cov : 0, lg_cov > 0 ? lg_cov - 1 : 0,
                sep_out, lg_cov >= 0 ? lg_t1 - lg_tb : 0.0, lg_t1 - lg_t0) ;
            fclose(f) ;
        }
    }

    return passes ;
}

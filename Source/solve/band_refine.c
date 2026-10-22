#include <stdlib.h>
#include <string.h>
#include "core.h"
#include "solve_internal.h"
#include "utils.h"
#include "vs.h"


// Two-sided band refinement of a vertex separator.
//
// The band is the separator S plus up to max_layers BFS layers into each
// shore. The rest of X is contracted into the source and the rest of Y into
// the sink. Each band vertex v is split into v_in -> v_out with capacity c(v),
// and graph edges inside the band get infinite capacity. A minimum s-t cut is
// then a minimum weight separator inside the band. Unlike a cover step, it can
// move into both shores at once, so it is not limited by the slack of one
// shore.
//
// The cut ignores the bound u on the shores. Of the two extreme minimum cuts
// (closest to the source and closest to the sink), the feasible one with the
// lightest separator, then the better balance, is taken. If neither is
// feasible, the pass makes no move. A move is accepted if the separator gets
// lighter, or stays the same weight and the shores get more balanced.

#define BAND_DEG_FACTOR 16  // total degree of the band <= this * max vertices

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
    int      *level, *it, *q, *stk ;
} Band_Net ;


static void Net_Init(Band_Net *N, int nn, long long max_ne)
{
    N->nn    = nn ;
    N->ne    = 0 ;
    N->head  = malloc(nn * sizeof(int)) ;
    N->level = malloc(nn * sizeof(int)) ;
    N->it    = malloc(nn * sizeof(int)) ;
    N->q     = malloc(nn * sizeof(int)) ;
    N->stk   = malloc(nn * sizeof(int)) ;
    N->nxt   = malloc(max_ne * sizeof(int)) ;
    N->to    = malloc(max_ne * sizeof(int)) ;
    N->cap   = malloc(max_ne * sizeof(band_cap)) ;
    for(int i = 0 ; i < nn ; i++) N->head[i] = -1 ;
}


static void Net_Free(Band_Net *N)
{
    free(N->head) ;
    free(N->level) ;
    free(N->it) ;
    free(N->q) ;
    free(N->stk) ;
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


// Dinic's algorithm on the residual capacities in N->cap
static void Net_Max_Flow(Band_Net *N, int s, int t)
{
    int *level = N->level ;
    int *it    = N->it ;
    int *stk   = N->stk ;
    int *q     = N->q ;

    while(true)
    {
        /* === Level graph ================================================== */
        int head = 0, tail = 0 ;
        for(int i = 0 ; i < N->nn ; i++) level[i] = -1 ;
        level[s] = 0 ;
        q[tail++] = s ;
        while(head < tail)
        {
            int u = q[head++] ;
            for(int e = N->head[u] ; e != -1 ; e = N->nxt[e])
            {
                int v = N->to[e] ;
                if(N->cap[e] > 0 && level[v] < 0)
                {
                    level[v]  = level[u] + 1 ;
                    q[tail++] = v ;
                }
            }
        }
        if(level[t] < 0) break ;

        /* === Blocking flow: iterative DFS with current-arc pointers ======= */
        for(int i = 0 ; i < N->nn ; i++) it[i] = N->head[i] ;
        int top = 0, u = s ;
        while(true)
        {
            if(u == t)
            {
                int      kmin = 0 ;
                band_cap f    = N->cap[stk[0]] ;
                for(int k = 1 ; k < top ; k++)
                {
                    if(N->cap[stk[k]] < f)
                    {
                        f    = N->cap[stk[k]] ;
                        kmin = k ;
                    }
                }
                for(int k = 0 ; k < top ; k++)
                {
                    N->cap[stk[k]]     -= f ;
                    N->cap[stk[k] ^ 1] += f ;
                }

                // Retreat to the tail of the bottleneck arc
                top = kmin ;
                u   = N->to[stk[kmin] ^ 1] ;
                continue ;
            }

            int e = it[u] ;
            while(e != -1 && !(N->cap[e] > 0 && level[N->to[e]] == level[u]+1))
            {
                e = N->nxt[e] ;
            }
            it[u] = e ;

            if(e != -1)
            {
                stk[top++] = e ;
                u = N->to[e] ;
            }
            else
            {
                // Dead end: remove u from the level graph and back up
                if(u == s) break ;
                level[u] = -1 ;
                e = stk[--top] ;
                u = N->to[e ^ 1] ;
                it[u] = N->nxt[it[u]] ;
            }
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

    Net_Max_Flow(N, S_NODE, T_NODE) ;
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


// Refine the separator X, Y with at most max_passes band moves. The band has
// at most 2 * max_bip vertices and max_layers BFS layers on each side. If the
// separator alone is too big for the band, or no move is accepted,
// Cover_Refine is run instead.
// Returns the number of band passes
int Band_Refine
(
    State   *X,
    State   *Y,
    Problem *P,
    VSPINT   max_bip,
    VSPINT   max_passes,
    VSPINT   max_layers
)
{
    SET_FUNCTION_COUNTER ;
    if(max_bip <= 0) return 0 ;

    VSPINT   i ;
    VSPLONG  p ;
    VSPINT   n        = P->graph->n ;
    VSPINT  *c        = P->graph->c ;
    VSPINT  *Ai       = P->graph->Ai ;
    VSPLONG *Ap       = P->graph->Ap ;
    VSPINT  *x        = X->point->int_part ;
    VSPINT  *y        = Y->point->int_part ;
    VSPINT   tot      = P->graph->total_vtx_cost ;
    VSPINT   max_band = 2 * max_bip ;
    int      passes   = 0 ;
    bool     moved    = false ;

    VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;

    int nS = 0 ;
    for(i = 0 ; i < n ; i++) if(!x[i] && !y[i]) nS++ ;
    if(nS == 0 || 2 * nS >= max_band)
    {
        Cover_Refine(X, Y, P, max_bip) ;
        return 0 ;
    }

    Band B ;
    B.x     = x ;
    B.y     = y ;
    B.c     = c ;
    B.u     = P->u ;
    B.band  = malloc(max_band * sizeof(int)) ;
    B.tmp   = malloc(max_band * sizeof(int)) ;
    B.best  = malloc(max_band * sizeof(int)) ;
    B.zs    = malloc(2 * max_band + 2) ;
    B.zt    = malloc(2 * max_band + 2) ;
    int *Slist = malloc(max_band * sizeof(int)) ;
    int *loc   = malloc(n * sizeof(int)) ;   // band index, -1 if not in band
    for(i = 0 ; i < n ; i++) loc[i] = -1 ;

    nS = 0 ;
    for(i = 0 ; i < n ; i++) if(!x[i] && !y[i]) Slist[nS++] = i ;

    while(passes < max_passes)
    {
        passes++ ;

        /* === Band: S plus max_layers BFS layers into each shore =========== */
        // Each shore gets at most half of the vertices left after S, and the
        // total degree of the band is capped at BAND_DEG_FACTOR * max_band
        int       nb       = 0, cnt[2] = {0, 0}, k ;
        int       cap_side = (max_band - nS) / 2 ;
        long long deg_sum  = 0 ;
        long long deg_cap  = (long long) BAND_DEG_FACTOR * max_band ;

        for(k = 0 ; k < nS ; k++)
        {
            loc[Slist[k]] = nb ;
            B.band[nb++]  = Slist[k] ;
            deg_sum      += Ap[Slist[k]+1] - Ap[Slist[k]] ;
        }
        if(deg_sum > deg_cap)
        {
            for(k = 0 ; k < nb ; k++) loc[B.band[k]] = -1 ;
            break ;
        }

        int lo = 0, hi = nb ;
        for(int layer = 0 ; layer < max_layers ; layer++)
        {
            for(k = lo ; k < hi ; k++)
            {
                int v = B.band[k] ;
                for(p = Ap[v] ; p < Ap[v+1] ; p++)
                {
                    int w = Ai[p] ;
                    if(loc[w] != -1) continue ;

                    int sd = SIDE(x, y, w) ;
                    if(sd == 2 || cnt[sd] >= cap_side) continue ;

                    long long dw = Ap[w+1] - Ap[w] ;
                    if(deg_sum + dw > deg_cap) continue ;

                    loc[w]       = nb ;
                    B.band[nb++] = w ;
                    cnt[sd]++ ;
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

        Net_Free(&B.N) ;
        for(k = 0 ; k < nb ; k++) loc[B.band[k]] = -1 ;

        if(!B.have_best) break ;

        /* === Accept ======================================================== */
        Band_Apply(&B, X, Y, P) ;
        moved = true ;
#ifdef DEBUG
        Band_Check_Separator(X, Y, P) ;
#endif

        // The new separator lies in the band
        nS = 0 ;
        for(k = 0 ; k < nb ; k++) if(B.best[k] == 2) Slist[nS++] = B.band[k] ;
        if(nS == 0 || 2 * nS >= max_band) break ;
    }

    free(B.band) ;
    free(B.tmp) ;
    free(B.best) ;
    free(B.zs) ;
    free(B.zt) ;
    free(Slist) ;
    free(loc) ;

    if(!moved) Cover_Refine(X, Y, P, max_bip) ;

    return passes ;
}

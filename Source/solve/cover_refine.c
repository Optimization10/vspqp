#include <stdlib.h>
#include "core.h"
#include "solve_internal.h"
#include "utils.h"
#include "vs.h"


// Number of vertices in the bipartite graph Cover_Step would build for shore
// `to`: separator vertices with a neighbor in the other shore, plus the
// other-shore vertices adjacent to the separator. This matches what
// CSC_Extract_Bipartite_Graph extracts. Counting stops once the count reaches
// limit, so the return value is only exact when it's below limit.
// mark is workspace of size n.
static VSPINT Bip_Size
(
    const int *where,
    int        to,
    int       *mark,
    VSPINT     limit,
    Problem   *P
)
{
    VSPINT   i, j ;
    VSPLONG  p ;
    VSPINT   nbip  = 0 ;
    VSPINT   n     = P->graph->n ;
    VSPINT  *Ai    = P->graph->Ai ;
    VSPLONG *Ap    = P->graph->Ap ;
    int      other = 1 - to ;

    for(i = 0 ; i < n ; i++) mark[i] = 0 ;

    for(i = 0 ; i < n && nbip < limit ; i++)
    {
        if(where[i] != 2) continue ;

        bool in_bip = false ;
        for(p = Ap[i] ; p < Ap[i+1] ; p++)
        {
            j = Ai[p] ;
            if(where[j] != other) continue ;

            in_bip = true ;
            if(!mark[j])
            {
                mark[j] = 1 ;
                nbip++ ;
            }
        }
        if(in_bip) nbip++ ;
    }

    return nbip ;
}


// One cover step toward shore `to` (0 = X, 1 = Y).
//
// where: current partition (0 = X, 1 = Y, 2 = separator), not modified
// lab:   output, the new partition in the same encoding
// Returns the new separator weight, or EMPTY if shore `to` would exceed u.
static VSPINT Cover_Step
(
    const int *where,
    int        to,
    int       *lab,
    Problem   *P
)
{
    VSPINT  i ;
    VSPINT  n     = P->graph->n ;
    VSPINT *c     = P->graph->c ;
    VSPINT  wto   = 0 ;
    VSPINT  wsep  = 0 ;
    int     other = 1 - to ;

    // vs_separator_min_weight covers the edges between label-0 and label-1
    // vertices and moves the cover to label 2. S goes on the flow's source
    // side (label 0) so that ties keep the current separator. Shore `to` gets
    // label 2 and is ignored.
    for(i = 0 ; i < n ; i++)
    {
        if     (where[i] == 2)     lab[i] = 0 ;
        else if(where[i] == other) lab[i] = 1 ;
        else                       lab[i] = 2 ;
    }

    vs_separator_min_weight(lab, P->graph->Ai, P->graph->Ap, n, c) ;

    // Decode back to the 0/1/2 partition encoding
    for(i = 0 ; i < n ; i++)
    {
        if     (where[i] == to) lab[i] = to ;    // shore `to` is unchanged
        else if(lab[i] == 0)    lab[i] = to ;    // S vertex not in the cover
        else if(lab[i] == 1)    lab[i] = other ; // `other` vertex not in cover
        // else lab[i] == 2: in the cover, i.e. the new separator

        if     (lab[i] == to) wto  += c[i] ;
        else if(lab[i] == 2)  wsep += c[i] ;
    }

    return (wto <= P->u) ? wsep : EMPTY ;
}


// Apply partition lab to X and Y. Only vertices whose side changed are
// touched, and Az is updated column by column
static void Apply_Partition(const int *lab, State *X, State *Y, Problem *P)
{
    VSPINT  i, s ;
    VSPINT  n         = P->graph->n ;
    VSPINT *c         = P->graph->c ;
    State  *states[2] = {X, Y} ;

    for(i = 0 ; i < n ; i++)
    {
        for(s = 0 ; s < 2 ; s++)
        {
            VSPINT *z    = states[s]->point->int_part ;
            VSPINT  want = (lab[i] == s) ;
            if(z[i] == want) continue ;

            z[i] = want ;
            Az_pm_col_of_A(states[s]->Az->int_part, i, want ? PLUS : MINUS, P) ;
            states[s]->weight += want ? c[i] : -c[i] ;
        }
    }
}


#ifdef DEBUG
// Check that X, Y form a true separator: no vertex in both shores and no edge
// between the shores. check_states covers the weights, Az, and the bound u.
// Returns the separator cost, recomputed from scratch.
static VSPINT Check_Separator(State *X, State *Y, Problem *P)
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

    // Error logs and continues, so report each kind of violation once
    if(nboth > 0)  ERR("Vertex in both shores") ;
    if(ncross > 0) ERR("Edge between the shores") ;
    if(wsep != P->graph->total_vtx_cost - X->weight - Y->weight)
    {
        ERR("Separator weight disagrees with the shore weights") ;
    }

    return wsep ;
}
#endif


// Repeat cover steps on both shores while the separator weight decreases.
// A direction is skipped when its bipartite graph has max_bip or more
// vertices
int Cover_Refine(State *X, State *Y, Problem *P, VSPINT max_bip)
{
    SET_FUNCTION_COUNTER ;
    if(max_bip <= 0) return 0 ;

    VSPINT  i, s, best ;
    VSPINT  n      = P->graph->n ;
    int     to, best_to ;
    int     passes = 0 ;
    int    *where  = malloc(n * sizeof(*where)) ;
    int    *lab[2] ;
    lab[0] = malloc(n * sizeof(int)) ;
    lab[1] = malloc(n * sizeof(int)) ;

    VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;
#ifdef DEBUG
    VSPINT wsep_prev = Check_Separator(X, Y, P) ;
#endif

    while(true)
    {
        passes++ ;

        for(i = 0 ; i < n ; i++)
        {
            if     (X->point->int_part[i]) where[i] = 0 ;
            else if(Y->point->int_part[i]) where[i] = 1 ;
            else                           where[i] = 2 ;
        }

        best    = P->graph->total_vtx_cost - X->weight - Y->weight ;
        best_to = EMPTY ;
        for(to = 0 ; to < 2 ; to++)
        {
            // lab[to] is scratch here; Cover_Step overwrites all of it
            if(Bip_Size(where, to, lab[to], max_bip, P) >= max_bip) continue ;

            s = Cover_Step(where, to, lab[to], P) ;
            if(s != EMPTY && s < best) { best = s ; best_to = to ; }
        }

        if(best_to == EMPTY) break ;
        Apply_Partition(lab[best_to], X, Y, P) ;

#ifdef DEBUG
        // An accepted step must leave a separator of the weight Cover_Step
        // reported, strictly below the previous one
        VSPINT wsep = Check_Separator(X, Y, P) ;
        if(wsep != best)      ERR("Separator weight differs from Cover_Step") ;
        if(wsep >= wsep_prev) ERR("Separator weight did not decrease") ;
        wsep_prev = wsep ;
#endif
    }

    free(where) ;
    free(lab[0]) ;
    free(lab[1]) ;

    return passes ;
}

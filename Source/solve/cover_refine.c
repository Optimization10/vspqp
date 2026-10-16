#include <stdlib.h>
#include "core.h"
#include "solve_internal.h"
#include "utils.h"
#include "vs.h"

/* ========================================================================== */
// Separator refinement by min-weight bipartite vertex cover
// (Liu 1989; Ashcraft & Liu 1998).
//
// Let S be the separator. For a shore `to` and the opposite shore `other`,
// take a min-weight vertex cover C of the bipartite graph between S and the
// vertices of `other` adjacent to S. Moving S \ C into `to` and C ∩ other into
// the separator gives a valid separator of weight c(C) <= c(S), since S is
// itself a cover.
/* ========================================================================== */


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
// touched, and Az is updated column by column, as in FM's Where_To_States.
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


// Repeat cover steps on both shores while the separator weight decreases.
// A direction is skipped when its bipartite graph has max_bip or more
// vertices; max_bip <= 0 disables the refinement.
//
// Requires X, Y to form a true separator with no fractional entries, e.g.
// right after FM_Metis. Terminates because each accepted step strictly
// decreases an integer separator weight.
void Cover_Refine(State *X, State *Y, Problem *P, VSPINT max_bip)
{
    if(max_bip <= 0) return ;

    VSPINT  i, s, best ;
    VSPINT  n     = P->graph->n ;
    int     to, best_to ;
    int    *where = malloc(n * sizeof(*where)) ;
    int    *lab[2] ;
    lab[0] = malloc(n * sizeof(int)) ;
    lab[1] = malloc(n * sizeof(int)) ;

    VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;

    while(true)
    {
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
    }

#ifdef DEBUG
    check_states(X, Y, P) ;
#endif

    free(where) ;
    free(lab[0]) ;
    free(lab[1]) ;
}

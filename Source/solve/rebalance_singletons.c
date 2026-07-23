#include "core.h"
#include "solve_internal.h"
#include "data_structures.h"
#include "utils.h"

#define PRINT 1

/* ========================================================================== */
// Whole-singleton insert/remove. Heap argument may be NULL (no heap insert).
/* ========================================================================== */

static void Insert_Singleton_Whole
(
    VSPINT       i,     // the singleton
    VSPINT       cval,  // cost of vertex i
    State       *state, // state to insert i into
    min_heap_v2 *h      // heap to insert i into (NULL to skip)
)
{
    state->point->int_part[i] = 1 ;
    state->Az->int_part[i]++ ;
    state->weight += cval ;
    if (h != NULL) MinHeap_Insert_v2(h, i) ;
}

// i has been removed from its heap outside this function
static void Remove_Singleton_Whole
(
    VSPINT       i,     // the singleton
    VSPINT       cval,  // cost of vertex i
    State       *state  // state to remove i from
)
{
    state->point->int_part[i] = 0 ;
    state->Az->int_part[i]-- ;
    state->weight -= cval ;
}


static bool Move_Singleton
(
    int to, State *states[2], min_heap_v2 *heaps[2], VSPINT u
)
{
    VSPINT cval, i, wto, woth, old_imb, new_imb ;
    VSPINT other = 1 - to ;

    // Check for empty heap
    MinHeap_See_Top_Key_v2(heaps[other], &cval) ;
    if(cval == EMPTY) return false ;

    // Check for infeasibility
    if(states[to]->weight + cval > u) return false ;

    // Check if new move improces the imbalance
    wto     = states[to]->weight;
    woth    = states[1 - to]->weight;    // the other shore
    old_imb = ABS_VAL(wto - woth) ;
    new_imb = ABS_VAL((wto + cval) - (woth - cval)) ;
    if (new_imb >= old_imb) return false ; // don't take a non-improving move

    // Move singleton from other to to
    MinHeap_Get_Top_Val_v2(heaps[other], &i) ;
    Remove_Singleton_Whole(i, cval, states[other]) ;
    Insert_Singleton_Whole(i, cval, states[to], heaps[to]) ;

    return true ;
}

static void Rebalance
(
    State *states[2], min_heap_v2 *heaps[2], Problem *P
)
{
    VSPINT to, moved ;
    bool stop = false ;
    while (!stop)
    {
        // Destination is the lighter state
        to    = (states[0]->weight < states[1]->weight) ? 0 : 1 ;
        moved = Move_Singleton(to, states, heaps, P->u) ;

        // Stop if there was no move, or if balance is nearly perfect
        if (!moved || imbalance(states[0], states[1]) <= 1) stop = true ;

        if(moved)
        {
#if PRINT
    char target = (to == 0) ? 'x' : 'y' ;
    printf("    Moved vertex. Imbal: %i\n", imbalance(states[0], states[1])) ;
#endif
        }
    }
}

// Round up or down a fractional entry, whichever way imrpoves the objective
// If objective doesn't change, round down
static void Round_Fractional_Entry(State *Z, Az *Aw, Problem *P)
{
    VSPINT u   = P->u ;
    VSPINT i   = Z->point->ind_frac ;
    VSPINT wgt = Z->weight ;

    if(i == EMPTY) return ;

    VSPINT ci = P->graph->c[i] ;

    VSPFLOAT g, a ;
    VSPFLOAT f = Z->point->val_frac ;
    VSPFLOAT gamma = P->graph->max_vtx_cost ;

    // Compute the partial derivative w.r.t. the fractional index
    g = ci - gamma * Azi(i, Aw) ;

    if(g > ZERO_PLUS)
    {
        // Round up

        // a is the increase in zi needed to hit the upper balance constraint
        a = (u - wgt) / ( (double) ci ) ;

        // we want to increase zi by the minimum of {a, 1-f}

        if(1 - f < a + ZERO_PLUS)
        {
            // The entire item fits
            Z->point->ind_frac = EMPTY ;
            Z->Az->ind_frac_z  = EMPTY ;

            Z->point->int_part[i] = 1 ;
            Z->weight += lround( (1-f) * ci ) ;

            Az_pm_col_of_A(Z->Az->int_part, i, PLUS, P) ;
        }
        else
        {
            Z->point->val_frac = f+ a ;
            Z->Az->val_frac_z  = f+ a ;
            Z->weight = u ;
        }
    }
    else
    {
        // Round down
        Z->point->ind_frac = EMPTY ;
        Z->Az->ind_frac_z  = EMPTY ;
        Z->weight -= lround(f * ci) ;
    }
}

// Assume i doesn't fit fully
static void Insert_Singleton_Fraction
(
    VSPINT   i,
    State   *state,
    Problem *P
)
{
    VSPINT one ;
    VSPINT u = P->u ;
    state->point->ind_frac = i ;
    state->point->val_frac = (u - state->weight) / (VSPFLOAT) (P->graph->c[i]) ;

    state->Az->ind_frac_z = i ;
    state->Az->val_frac_z = state->point->val_frac ;

    one = ++(state->Az->one) ;
    state->Az->col_frac_z[i] = one ;

    state->weight = u ;
}

// If a state already has a fractional entry, this function doesn't touch it
static void Use_Separator_Heap
(
    State       *states[2],
    min_heap_v2 *heap,
    Problem     *P
)
{
    if(heap == NULL) return ;
    int to, other ;
    VSPINT i, ci ;
    VSPINT *c = P->graph->c ;
    VSPINT u  = P->u ;

    // Stop when the heap is empty ;
    while(true)
    {
        // Pop off member of heap. Quit if heap is empty
        MinHeap_Get_Top_Val_v2(heap, &i) ;
        if(i == EMPTY) return ;

        ci = c[i] ;

        // Determine smaller shore
        to = (states[0]->weight < states[1]->weight) ? 0 : 1 ;
        other = 1 - to ;

        // Add i to smaller shore whole if possible, or a fraction if that
        // shore doesn't already have a fraction
        if(states[to]->weight + ci <= u)
        {
            Insert_Singleton_Whole(i, ci, states[to], NULL) ;
        }
        else if(states[to]->weight < u && states[to]->point->ind_frac == EMPTY)
        {
            Insert_Singleton_Fraction(i, states[to], P) ;
        }
        else if(states[other]->weight < u
                && states[other]->point->ind_frac == EMPTY)
        {
            // Try putting in the heavier shore as a last ditch attempt
            Insert_Singleton_Fraction(i, states[other], P) ;
        }
    } 
}

void Rebalance_Singletons(State *X, State *Y, Problem *P)
{
    SET_FUNCTION_COUNTER ;
#ifdef DEBUG
    check_states(X, Y, P) ;
#endif

#if PRINT
    printf("Problem: %s  Level: %i  n: %i  edges: %i\n", G_Current_Problem_Name, P->level, P->graph->n, P->graph->nz/2) ;
#endif

    /* ====================================================================== */
    // Quit early if partition is already near-perfectly balanced
    /* ====================================================================== */
    if (imbalance(X, Y) <= 1) return ;

    /* ====================================================================== */
    // Declare variables
    /* ====================================================================== */
    VSPINT i, to ;
    VSPINT   u  = P->u ;
    VSPINT  *c  = P->graph->c ;
    VSPINT   n  = P->graph->n ;
    VSPLONG *Ap = P->graph->Ap ;
    State *states[2] = {X, Y} ;
    VSPINT *z[2] = {X->point->int_part, Y->point->int_part} ;
    VSPINT idx_frac[2] = {X->point->ind_frac, Y->point->ind_frac} ;
    bool moved ;

    /* ====================================================================== */
    // Build the two heaps, one for x and one for y
    /* ====================================================================== */
    min_heap_v2 *heaps[2] ;
    min_heap_v2 *sep_heap = NULL ;
    heaps[0] = MinHeap_Init_v2(c, n) ;
    heaps[1] = MinHeap_Init_v2(c, n) ;

    // Let i be a singleton. If i is in X only, it gets added to the x heap,
    // and if i is in Y only, it gets added to the Y heap. If i is
    // At the end of this loop:
    // 1) All singletons fully in a shore are in the heap for that shore
    // 2) No singletons are in both shores (not even fractionally. Fractional
    //    entries of both variables are skipped)
    // 3) Any singleton in the separator is in the heap corresponding to the
    //    separator
    // The code many have modified the original point in the following way:
    // It may have moved singletons from the separator to the shores,
    // and it moved any singletons in the intersection to just one shore
    for (i = 0 ; i < n ; i++)
    {
        // Skip non-singletons
        if (Ap[i+1] - Ap[i] > 0) continue ;

        // Skip fractional entries
        if (i == idx_frac[0] || i == idx_frac[1])
        {

#if PRINT
    printf("\n Problem %s .Fractional Singleton", G_Current_Problem_Name) ;
#endif
            continue ;
        }
        

        // Check if singleton is in intersection. Really shouldn't happen
        // (it incurs penalty gamma >= c[i])
        if (z[0][i] + z[1][i] >= 2)
        {
            printf("WARNING: Singleton in intersection\n") ;
            to = (states[0]->weight > states[1]->weight) ? 0 : 1 ;
            Remove_Singleton_Whole(i, c[i], states[to]) ;
        }

        // Separator singleton. Very unlikely. Try inserting it whole into
        // the smaller shore. If not possible, then place it into a separator
        // heap
        if (z[0][i] + z[1][i] == 0)
        {
#if PRINT
    printf("\n Problem %s . Singleton in separator", G_Current_Problem_Name) ;
#endif
            to = (states[0]->weight < states[1]->weight) ? 0 : 1 ;
            if (states[to]->weight  + c[i] <= u)
            {
                // Insert singleton into shore and into heap
                Insert_Singleton_Whole(i, c[i], states[to], NULL) ;
                heaps[to]->heap[ ++(heaps[to]->numNodes) ] = i ;
            }
            else
            {
                if(sep_heap == NULL) sep_heap = MinHeap_Init_v2(c, n) ;
                MinHeap_Insert_v2(sep_heap, i) ;
            }
            continue ;
        }

        // Append i to the heap of the shore that contains it
        to = (states[0]->point->int_part[i]) ? 0 : 1 ;
        heaps[to]->heap[ ++(heaps[to]->numNodes) ] = i ;
    }

    /* ====================================================================== */
    // Quit early if no singleton in the heavier shore
    /* ====================================================================== */
    to = (states[0]->weight < states[1]->weight) ? 1 : 0 ;
    if (heaps[to]->numNodes == 0) goto CLEANUP ;

    /* ====================================================================== */
    // Heapify the arrays
    /* ====================================================================== */
    MinHeap_Heapify_v2(heaps[0]) ;
    MinHeap_Heapify_v2(heaps[1]) ;

    /* ====================================================================== */
    // Perform iterations
    /* ====================================================================== */
    Rebalance(states, heaps, P) ;

#ifdef DEBUG
        check_states(X, Y, P) ;
#endif

    /* ====================================================================== */
    // Cleanup
    /* ====================================================================== */
CLEANUP:
    // Use up the separator heap
    Use_Separator_Heap(states, sep_heap, P) ;

    Round_Fractional_Entry(X, Y->Az, P) ;
    Round_Fractional_Entry(Y, X->Az, P) ;
    for(i = 0 ; i < 2 ; i++) MinHeap_Free_v2(heaps[i]) ;
    if(sep_heap != NULL) MinHeap_Free_v2(sep_heap) ;

#ifdef DEBUG
    check_states(X, Y, P) ;
#endif
    return ;
}
#include "solve_internal.h"
#include "data_structures.h"
#include "utils.h"

static void Insert_Singleton_Whole
(
    VSPINT       i,     // the singleton
    VSPINT       cval,  // cost of vertex i
    State       *state, // state to insert i into
    min_heap_v2 *h      // heap to insert i into
)
{
    state->point->int_part[i] = 1 ;
    state->Az->int_part[i]++ ;
    state->weight += cval ;
    MinHeap_Insert_v2(h, i) ;
}

// i has been removed from its heap outside this function
static void Remove_Singleton_Whole
(
    VSPINT       i,     // the singleton
    VSPINT       cval,  // cost of vertex i
    State       *state // state to remove i from
)
{
    state->point->int_part[i] = 0 ;
    state->Az->int_part[i]-- ;
    state->weight -= cval ;
}

// Assumed: the singleton does not fit wholly into the shore
static void Insert_Singleton_Fractional
(
    VSPINT       i,     // the singleton
    VSPFLOAT     cval,  // cost of vertex i
    State       *state, // state to insert i into
    VSPINT       u
)
{
    VSPFLOAT val_frac = (u - state->weight) / cval ;

    // Update point
    state->point->int_part[i] = 0 ;
    state->point->ind_frac = i ;
    state->point->val_frac = val_frac ;

    // Update Az
    state->Az->ind_frac_z  = i ;
    state->Az->val_frac_z  = val_frac ;
    state->Az->one++ ;
    state->Az->col_frac_z[ i ] = state->Az->one ;

    // Update cost
    state->weight = u ;
}

// Assume that i is a singleton and the fractional entry in the state
static void Remove_Singleton_Fractional
(
    VSPINT       i,     // the singleton
    VSPFLOAT     cval,  // cost of vertex i
    State       *state // state to remove i from
)
{
    VSPFLOAT val_frac = state->point->val_frac ;

    // Update point
    state->point->ind_frac = EMPTY ;
    state->point->val_frac = 0 ;

    // Update Az
    state->Az->ind_frac_z  = EMPTY ;
    state->Az->val_frac_z  = 0 ;

    // Update cost
    state->weight -= lround( cval * val_frac ) ;
}

// Attempt to move the root of the source heap (source = 0, 1, or 2) to the
// to state (to = 0 or 1)
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


// An iteration of the algorithm consts of moving a singleton
// from the separator or heavier shore (separator preferred)
// The iterations stop when one of the following stopping conditions are met:
// Stopping conditions:
// 1) Difference in weight is at most 1
// 2) No singletons in separator or heavier variable
// 3) Each singleton has moved N number of times (N = something like 2)
// 4) Moving any singleton from separator or heavier shore to lighter one is
//    infeasible
//
// Implementation detail: the third criterion is implicit. After a singleton
// has moved its maximum number of times, it is simply not added to the heap
// Make sure that if there is a fractional entry, then the cost of the shore
// is u, as this is something expected by the rest of the code
void Rebalance_Singletons(State *X, State *Y, Problem* P)
{
    SET_FUNCTION_COUNTER ;
#ifdef DEBUG

    check_states(X, Y, P) ;
#endif
    /* ====================================================================== */
    // Variable declarations
    /* ====================================================================== */
    VSPINT i, j, to, other, cval ;
    VSPINT  u          = P->u ;
    VSPINT *c          = P->graph->c ;
    VSPINT  n          = P->graph->n ;
    VSPINT *z[2]       = {X->point->int_part, Y->point->int_part} ;
    VSPINT idx_frac[2] = {X->point->ind_frac, Y->point->ind_frac} ;

    VSPLONG *Ap = P->graph->Ap ;
    State *states[2] = {X, Y} ;  
    bool stop, moved ; 

    /* ====================================================================== */
    // Quit early if partition is already near-perfectly balanced
    /* ====================================================================== */
    if( imbalance(X, Y) <= 1 ) return ;

    /* ====================================================================== */
    // Build the 2 min heaps. 0 is for x, 1 for y
    /* ====================================================================== */
    min_heap_v2 *heaps[2] ;
    for(i = 0 ; i < 2 ; i++) heaps[i] = MinHeap_Init_v2(c, n) ;

    // Insert singletons into the x heap, y heap
    for(i = 0 ; i < n ; i++)
    {
        // Skip non-singletons
        if(Ap[i+1] - Ap[i] > 0) continue ;

        // Check if singleton is in intersection. Really shouldn't happen
        if(z[0][i] + z[1][i] >= 2)
        {
            printf("WARNING: Singleton in intersection\n") ;

            // Remove the singleton from the bigger shore
            to = (states[0]->weight > states[1]->weight) ? 0 : 1 ;
            Remove_Singleton_Whole(i, c[i], states[to]) ;
        }

        // Skip singletons that have a fractional value
        if(i == idx_frac[0] || i == idx_frac[1]) continue ;

        // Place i into the appropriate heap
        to = 2 ;
        for(j = 0 ; j < 2 ; j++)
        {
            if( states[j]->point->int_part[i] ) to = j ;
        }

        // If by some miracle to is in the separator, just leave it alone
        if(to != 2)
        {
            MinHeap_Insert_v2(heaps[to], i) ;
        }
    }

    /* ====================================================================== */
    // Quit early if no singleton in the heavier shore
    /* ====================================================================== */
    // to is the bigger shore
    to = (states[0]->weight < states[1]->weight ? 1 : 0) ;
    if(heaps[to]->numNodes == 0) goto FRACTIONS ;

    /* ====================================================================== */
    // Perform iterations
    /* ====================================================================== */
    stop = false ;

    while( !stop )
    {
        // Identify the smallest state
        to    = (states[0]->weight < states[1]->weight) ? 0 : 1 ;

        // Shift singleton from larger to smaller state
        moved = Move_Singleton(to, states, heaps, u) ;

        // Stop if there was no move, or if balance is nearly perfect
        if(!moved || imbalance(states[0], states[1]) <= 1) stop = true ;

#ifdef DEBUG
    check_states(X, Y, P) ;
#endif
    }
    
FRACTIONS:
    /* ====================================================================== */
    // Deal with fractional entries
    // If the bigger shore has a fractional entry which is a singleton,
    // consider transferring it to the smaller shore
    /* ====================================================================== */
    to    = (states[0]->weight < states[1]->weight) ? 0 : 1 ;
    other = 1 - to ;
    if(idx_frac[other] != EMPTY)
    {
        i = idx_frac[other] ;
        if(Ap[i+1] - Ap[i] == 0)
        {
            // Check if transferring i from the heavier state and placing it
            // fully in the lighter state is feasible and improves imbalance
            cval = c[i] ;
            VSPFLOAT val_frac    = states[other]->point->val_frac ;
            VSPINT   frac_wgt    = (VSPINT) lround(cval * val_frac) ;
            VSPINT   new_wgt_to    = states[to]->weight + c[i] ;
            VSPINT   new_wgt_other = states[other]->weight - frac_wgt ;

            VSPINT old_imb = ABS_VAL(states[to]->weight - states[other]->weight) ;
            VSPINT new_imb = ABS_VAL(new_wgt_to - new_wgt_other) ;

            if(new_wgt_to <= u && new_imb < old_imb
               && states[to]->point->int_part[i] == 0
               && states[to]->point->ind_frac != i )
            {
                // Place i fully into the lighter state
                states[to]->point->int_part[i] = 1 ;
                states[to]->Az->int_part[i]++ ;
                states[to]->weight += cval ;

                // Remove the fractional entry from the heavier state
                Remove_Singleton_Fractional(i, cval, states[other]) ;
            }
        }
    }


    /* ====================================================================== */
    // Free memory
    /* ====================================================================== */
    for(i = 0 ; i < 2 ; i++) MinHeap_Free_v2(heaps[i]) ;

#ifdef DEBUG
    check_states(X, Y, P) ;
#endif
}

#include "core.h"
#include <stdlib.h>
#include "hubs.h"

// Return a list of vertices whose degree / n is > t
void Get_Hubs
(
    VSPINT      n,
    VSPLONG    *Ap,
    VSPFLOAT    t,
    VSPINT    **hubs,    // Return
    VSPINT     *num_hubs // Return
)
{
    VSPINT i ;
    VSPFLOAT n_reciprocal = 1 / (VSPFLOAT) n ;

    *hubs     = malloc(n * sizeof(**hubs)) ;
    *num_hubs = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        if((Ap[i+1] - Ap[i]) * n_reciprocal > t)
        {
            (*hubs)[ (*num_hubs)++ ] = i ;
        }
    }
}

void Remove_Hubs
(
    VSPINT    n_input,
    VSPINT   *Ai_input,
    VSPLONG  *Ap_input,
    VSPINT   *hubs,
    VSPINT    num_hubs,
    VSPINT   *n,    // output
    VSPINT  **Ai,   // output
    VSPLONG **Ap    // output
)
{
    VSPINT new_n = n_input - num_hubs ;
    *n = new_n ;

    // Build a map from old vertex index to new vertex index.  Hubs are
    // flagged with EMPTY. Every surviving vertex is renumbered 0, 1, 2, ...
    VSPINT *old_to_new = (VSPINT *) malloc (n_input * sizeof (VSPINT)) ;
    for (VSPINT v = 0 ; v < n_input ; v++) old_to_new [v] = 0 ;
    for (VSPINT k = 0 ; k < num_hubs ; k++) old_to_new [hubs [k]] = EMPTY ;
    VSPINT count = 0 ;
    for (VSPINT v = 0 ; v < n_input ; v++)
    {
        if (old_to_new [v] != EMPTY) old_to_new [v] = count++ ;
    }

    // Allocate the new column-pointer array.
    VSPLONG *Ap_new = (VSPLONG *) malloc ((new_n + 1) * sizeof (VSPLONG)) ;

    // First pass: count the surviving entries in each new column.  Counts
    // are stored shifted by one so the prefix sum below produces Ap_new.
    for (VSPINT j = 0 ; j <= new_n ; j++) Ap_new [j] = 0 ;
    for (VSPINT j = 0 ; j < n_input ; j++)
    {
        if (old_to_new [j] == EMPTY) continue ;            // skip hub columns
        VSPLONG nz = 0 ;
        for (VSPLONG p = Ap_input [j] ; p < Ap_input [j+1] ; p++)
        {
            if (old_to_new [Ai_input [p]] != EMPTY) nz++ ;  // keep non-hub rows
        }
        Ap_new [old_to_new [j] + 1] = nz ;
    }

    // Prefix sum: turn per-column counts into column pointers.
    for (VSPINT j = 0 ; j < new_n ; j++) Ap_new [j+1] += Ap_new [j] ;

    // Allocate the new row-index array.
    VSPLONG nnz = Ap_new [new_n] ;
    VSPINT *Ai_new = (VSPINT *) malloc (nnz * sizeof (VSPINT)) ;

    // Second pass: fill in the renumbered row indices.
    for (VSPINT j = 0 ; j < n_input ; j++)
    {
        if (old_to_new [j] == EMPTY) continue ;
        VSPLONG p_out = Ap_new [old_to_new [j]] ;
        for (VSPLONG p = Ap_input [j] ; p < Ap_input [j+1] ; p++)
        {
            VSPINT i = Ai_input [p] ;
            if (old_to_new [i] != EMPTY) Ai_new [p_out++] = old_to_new [i] ;
        }
    }

    free (old_to_new) ;

    *Ai = Ai_new ;
    *Ap = Ap_new ;
}


void Restore_Hubs
(
    State    *state,     // modified in place (its point only)
    VSPINT   *hubs,
    VSPINT    num_hubs
)
{
    Point *point = state->point ;

    VSPINT  n_small = point->n ;             // size of the reduced graph
    VSPINT  n_large = n_small + num_hubs ;   // size of the original graph

    VSPINT *int_small = point->int_part ;
    VSPINT *int_large = (VSPINT *) malloc (n_large * sizeof (VSPINT)) ;

    // Iterate through the original index space.  At each original vertex,
    // either it is the next hub (insert a 0 and advance the hub pointer) or it
    // is the next surviving vertex (copy its value from the reduced point).
    // Hubs are sorted in increasing order
    VSPINT new_frac = EMPTY ;   // fractional index in the original numbering
    VSPINT h = 0 ;              // pointer into hubs
    VSPINT s = 0 ;              // pointer into the reduced point

    // Each entry in the large point corresponds either to a hub, or to
    // an entry in the small point
    for (VSPINT v = 0 ; v < n_large ; v++)
    {
        if (h < num_hubs && hubs [h] == v)
        {
            int_large [v] = 0 ; // hub: outside the subset
            h++ ;
        }
        else
        {
            if (s == point->ind_frac) new_frac = v ;   // remap fractional index
            int_large [v] = int_small [s] ;
            s++ ;
        }
    }

    free (int_small) ;

    point->int_part = int_large ;
    point->ind_frac = new_frac ;   // EMPTY stays EMPTY (never matched above)
    point->n        = n_large ;
    // point->val_frac is unchanged.
    // state->weight and state->Az are left untouched.
}

void Place_Hubs_In_Shores
(
    State *X, State *Y,
    VSPINT *hubs, VSPINT num_hubs,
    VSPINT *Ai_input, VSPLONG *Ap_input,   // ORIGINAL graph
    VSPINT u                               // balance bound
)
{
    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;

    for(VSPINT k = 0 ; k < num_hubs ; k++)
    {
        VSPINT h = hubs[k] ;

        bool seen_x = false, seen_y = false ;
        for(VSPLONG p = Ap_input[h] ; p < Ap_input[h + 1] ; p++)
        {
            VSPINT w = Ai_input[p] ;
            if(x[w]) seen_x = true ;
            if(y[w]) seen_y = true ;
            if(seen_x && seen_y) break ; // touches both shores -> stuck
        }

        if(seen_x && seen_y) continue ; // stays in the separator

        if(seen_x) // X-side (and/or separator) only
        {
            if(X->weight + 1 <= u) { x[h] = 1 ; X->weight += 1 ; }
        }
        else if(seen_y) // Y-side (and/or separator) only
        {
            if(Y->weight + 1 <= u) { y[h] = 1 ; Y->weight += 1 ; }
        }
        else // no shore neighbors
        {
            bool x_ok = (X->weight + 1 <= u) ;
            bool y_ok = (Y->weight + 1 <= u) ;
            if(x_ok && (!y_ok || X->weight <= Y->weight))
            {
                x[h] = 1 ;
                X->weight += 1 ;
            }
            else if(y_ok)
            {
                y[h] = 1 ; Y->weight += 1 ;
            }
            // neither fits -> stays in the separator
        }
    }
}
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
    // flagged with -1; every surviving vertex is renumbered 0, 1, 2, ...
    // in its original relative order.
    VSPINT *old_to_new = (VSPINT *) malloc (n_input * sizeof (VSPINT)) ;
    for (VSPINT v = 0 ; v < n_input ; v++) old_to_new [v] = 0 ;
    for (VSPINT k = 0 ; k < num_hubs ; k++) old_to_new [hubs [k]] = -1 ;
    VSPINT count = 0 ;
    for (VSPINT v = 0 ; v < n_input ; v++)
    {
        if (old_to_new [v] != -1) old_to_new [v] = count++ ;
    }

    // Allocate the new column-pointer array.
    VSPLONG *Ap_new = (VSPLONG *) malloc ((new_n + 1) * sizeof (VSPLONG)) ;

    // First pass: count the surviving entries in each new column.  Counts
    // are stored shifted by one so the prefix sum below produces Ap_new.
    for (VSPINT j = 0 ; j <= new_n ; j++) Ap_new [j] = 0 ;
    for (VSPINT j = 0 ; j < n_input ; j++)
    {
        if (old_to_new [j] == -1) continue ;            // skip hub columns
        VSPLONG nz = 0 ;
        for (VSPLONG p = Ap_input [j] ; p < Ap_input [j+1] ; p++)
        {
            if (old_to_new [Ai_input [p]] != -1) nz++ ;  // keep non-hub rows
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
        if (old_to_new [j] == -1) continue ;
        VSPLONG p_out = Ap_new [old_to_new [j]] ;
        for (VSPLONG p = Ap_input [j] ; p < Ap_input [j+1] ; p++)
        {
            VSPINT i = Ai_input [p] ;
            if (old_to_new [i] != -1) Ai_new [p_out++] = old_to_new [i] ;
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

    // Walk the original index space.  At each original vertex, either it is
    // the next hub (insert a 0 and advance the hub pointer) or it is the next
    // surviving vertex (copy its value from the reduced point).  Hubs are
    // sorted in increasing order, so a single merge walk suffices.
    VSPINT new_frac = EMPTY ;   // fractional index in the original numbering
    VSPINT h = 0 ;              // pointer into hubs
    VSPINT s = 0 ;              // pointer into the reduced point
    for (VSPINT v = 0 ; v < n_large ; v++)
    {
        if (h < num_hubs && hubs [h] == v)
        {
            int_large [v] = 0 ;                        // hub: outside the subset
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
#include "core.h"
#include "data_structures.h"

#define FIRST_SMALLER(a, b, c, d) ((a)*(d) < (c)*(b))
#define AZ(i) (Az_int[ i ] + val_frac_z * (is_z_frac) * (one == col_A[ i ]) )


/******************************************************************************
VERSION 3
Sort by (Az)_i / c_i
*******************************************************************************/

/* On input:
The subtree of "heap" starting at node with index "p" is a heap,
except that it's root may be out of order.

On output: 
subtree of "heap" starting at node with index "p" will be a subheap.
*/
void min_heapify_v3
(
    VSPINT   *heap, /* size n + 1, sorted by cost array X */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*      c, /* size n. Denominator vector to sort by */
    VSPINT       p, /* start at node p in the heap */
    VSPINT   nheap  /* heap [1 ... n] is in use */
)
{
    VSPINT left, right, e, hleft, hright ;

    e = heap [p] ;
    VSPINT *Az_int = Az->int_part ;
    VSPINT ind_frac_z = Az->ind_frac_z ;
    VSPFLOAT val_frac_z = Az->val_frac_z ;
    VSPINT is_z_frac = (ind_frac_z != EMPTY) ;
    VSPINT *col_A = Az->col_frac_z ;
    VSPINT one = Az->one ;

    VSPFLOAT num_e = AZ(e) ;
    VSPINT denom_e = c[e] ;

    while ( TRUE )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            /* both left and right children of p are in the heap */
            hleft = heap [left] ;
            hright = heap [right] ;
            // Xleft = X [hleft] ;
            // Xright = X [hright] ;
            // if ((Xleft < Xright)) /* for min heap: (Xleft < Xright) */
            if ( FIRST_SMALLER( AZ(hleft), c[hleft], AZ(hright), c[hright]) )
            {
                /* traverse down the left heap */
                // if ((Xe > Xleft)) /* for min heap: (Xe > Xleft) */
                if ( FIRST_SMALLER(AZ(hleft), c[hleft], num_e, denom_e) )
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
            else
            {
                /* traverse down the right heap */
                // if ((Xe > Xright)) /* for min heap: (Xe > Xright) */
                if (FIRST_SMALLER(AZ(hright), c[hright], num_e, denom_e)) /* for min heap: (Xe > Xright) */
                {
                    /* swap e with the right child */
                    heap [p] = hright ;
                    p = right ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
        }
        else
        {
            /* only the left child of p is in the heap */
            if ( left <= nheap )
            {
                /* traverse down the left heap */
                hleft = heap [left] ;
                // Xleft = X [hleft] ;
                // if ((Xe > Xleft)) /* for min heap: (Xe > Xleft) */
                if ( FIRST_SMALLER(AZ(hleft), c[hleft], num_e, denom_e) )
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
            }
            /* element e, with value Xe, has found its final home */
            heap [p] = e ;
            return ;
        }
    }
}


/*** Version 3: Heapify by a the value of a ratio, Azi/ci */
void min_heap_build_v3
(
    VSPINT*    heap, /* on input, an unsorted set of indices */
    Az*      Az, /* size n. Numerator vector to sort by */
    VSPINT*       c, /* size n. Denominator vector to sort by */
    VSPINT        n  /* number of elements to build into the heap */
)
{
    VSPINT p  ;

    for (p = n/2  ; p >= 1  ; p--)
    {
        min_heapify_v3 (heap, Az, c, p, n)  ;
    }

}


/*
min_heap_add
On input:
1) array "heap" is a heap
2) "nheap" is number of nodes in heap

On output:
1) "heap" recieves "new" and is heapified 
2) "nheap" increases by 1
*/
void min_heap_add_v3 /* returns the new size of the heap */
(
    VSPINT *heap, /* size n + 1. containing indices into X */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT new, /* the new new */
    VSPINT *nheap /* number of elements in heap not counting new one */
)
{
    // VSPFLOAT num_old, Xnew ;
    VSPINT num_old, denom_old ;
    VSPINT i, parent, current ;

    VSPINT *Az_int = Az->int_part ;
    VSPINT ind_frac_z = Az->ind_frac_z ;
    VSPFLOAT val_frac_z = Az->val_frac_z ;
    VSPINT is_z_frac = (ind_frac_z != EMPTY) ;
    VSPINT *col_A = Az->col_frac_z ;
    VSPINT one = Az->one ;


    (*nheap)++ ;
    current = *nheap ;
    heap [current] = new ;
    // Xold = X [new] ;
    num_old = AZ(new) ;
    denom_old = c[new] ;

    while ( current > 1 )
    {
        parent = current/2 ;
        i = heap [parent] ;
        // Xnew = X [i] ;
        // if ((Xnew > Xold)) /* for min heap: ( Xnew > Xold ) */
        if ( FIRST_SMALLER(num_old, denom_old, AZ(i), c[i]) )
        {
            /* Xnew and Xold are out of order, so swap them */
            heap [parent] = new ;
            heap [current] = i ;
        }
        else return;
        current = parent ;
    }
    return;
}


/*
min_heap_delete

On input:
1) array "heap" is a heap
2) "nheap" is number of nodes in heap

On output:
1) array "heap" looses root and is heapified
2) "nheap", reduced by 1
*/
void min_heap_delete_v3 /* return new size of heap */
(
    VSPINT *heap, /* indices into X, 1..n on input.  size nheap+1 */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT *nheap /* number of items in heap */
)
{
    if (*nheap <= 1)
    {
        *nheap = 0;
        return;
    }

    /* move element from the end of the heap to the top */
    heap [1] = heap [*nheap] ;
    (*nheap)--;

    /* Heapify */
    min_heapify_v3 (heap, Az, c, 1, *nheap) ;
    return;
}
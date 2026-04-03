#include "core.h"
#include "shared_api.h"
#include "ds_api.h"



/******************************************************************************
VERSION 1
Heapify by a float array called X
*******************************************************************************/


/* On input:
The subtree of "heap" starting at node with index "p" is a heap,
except that it's root may be out of order.

On output: 
subtree of "heap" starting at node with index "p" will be a subheap.
*/
void min_heapify_v1
(
    VSPINT            *heap, /* size n + 1, sorted by cost array X */
    VSPFLOAT const       *X, /* size n. Cost vector to sort by */
    VSPINT                p, /* start at node p in the heap */
    VSPINT            nheap  /* heap [1 ... n] is in use */
)
{
    VSPFLOAT Xe, Xleft, Xright ;
    VSPINT left, right, e, hleft, hright ;

    e = heap [p] ;
    Xe = X [e] ;

    while ( TRUE )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            /* both left and right children of p are in the heap */
            hleft = heap [left] ;
            hright = heap [right] ;
            Xleft = X [hleft] ;
            Xright = X [hright] ;
            if ((Xleft < Xright)) /* for min heap: (Xleft < Xright) */
            {
                /* traverse down the left heap */
                if ((Xe > Xleft)) /* for min heap: (Xe > Xleft) */
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
                if ((Xe > Xright)) /* for min heap: (Xe > Xright) */
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
                Xleft = X [hleft] ;
                if ((Xe > Xleft)) /* for min heap: (Xe > Xleft) */
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

/*
min_heap_add
On input:
1) array "heap" is a heap
2) "nheap" is number of nodes in heap

On output:
1) "heap" recieves "leaf" and is heapified 
2) "nheap" increases by 1
*/
void min_heap_add_v1 /* returns the new size of the heap */
(
    VSPINT      *heap, /* size n + 1. containing indices into X */
    VSPFLOAT const *X, /* size n. Cost vector to sort by */
    VSPINT       leaf, /* the new leaf */
    VSPINT     *nheap /* number of elements in heap not counting new one */
)
{
    VSPFLOAT Xold, Xnew ;
    VSPINT i, newindex, old ;

    (*nheap)++ ;
    old = *nheap ;
    heap [old] = leaf ;
    Xold = X [leaf] ;
    while ( old > 1 )
    {
        newindex = old/2 ;
        i = heap [newindex] ;
        Xnew = X [i] ;
        if ((Xnew > Xold)) /* for min heap: ( Xnew > Xold ) */
        {
            /* Xnew and Xold are out of order, so swap them */
            heap [newindex] = leaf ;
            heap [old] = i ;
        }
        else return;
        old = newindex ;
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
void min_heap_delete_v1 /* return new size of heap */
(
    VSPINT           *heap, /* indices into X, 1..n on input.  size nheap+1 */
    VSPFLOAT const      *X, /* not modified */
    VSPINT          *nheap /* number of items in heap */
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
    min_heapify_v1 (heap, X, 1, *nheap) ;
    return;
}
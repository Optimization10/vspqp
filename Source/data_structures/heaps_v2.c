/******************************************************************************
VERSION 2 (struct-based, with lookup)

    min_heap_v2 / max_heap_v2 bundle:
        heap     : 1-indexed array of indices into c (owned by the struct)
        lookup   : lookup[i] = heap position of index i, or EMPTY (owned)
        c        : the key array (borrowed, not owned, not freed)
        n        : length of c (= capacity of the heap)
        numNodes : current number of elements

    Invariant: heap[lookup[i]] == i for every i in the heap;
               lookup[i] == EMPTY for every i not in the heap.
*******************************************************************************/

#include "core.h"
#include "data_structures.h"
#include "utils.h"
#include <stdlib.h>

/*----------------------------------------------------------------------------*/
// Max heap
/*----------------------------------------------------------------------------*/

// Allocates the heap object and its internal arrays.
// Returns NULL on allocation failure.
// c must have length n and must outlive the heap.
max_heap_v2 *MaxHeap_Init_v2
(
    VSPINT *c,
    VSPINT  n
)
{
    max_heap_v2 *h = malloc(sizeof(max_heap_v2)) ;
    if (h == NULL)
    {
        return NULL ;
    }

    h->heap   = malloc((size_t)(n + 1) * sizeof(VSPINT)) ;
    h->lookup = malloc((size_t) n      * sizeof(VSPINT)) ;
    if (h->heap == NULL || h->lookup == NULL)
    {
        free(h->heap) ;
        free(h->lookup) ;
        free(h) ;
        return NULL ;
    }

    h->c        = c ;
    h->n        = n ;
    h->numNodes = 0 ;

    for (VSPINT i = 0 ; i < n ; i++)
    {
        h->lookup[i] = EMPTY ;
    }
    return h ;
}

// Frees the internal arrays AND the heap object itself.
// Safe to call with NULL. The pointer is invalid after this call.
void MaxHeap_Free_v2
(
    max_heap_v2 *h
)
{
    if (h == NULL)
    {
        return ;
    }
    free(h->heap) ;
    free(h->lookup) ;
    free(h) ;
    return ;
}

// Empty the heap without freeing (O(n): re-EMPTYs the lookup array)
void MaxHeap_Reset_v2
(
    max_heap_v2 *h
)
{
    h->numNodes = 0 ;
    for (VSPINT i = 0 ; i < h->n ; i++)
    {
        h->lookup[i] = EMPTY ;
    }
    return ;
}

// O(1) membership test: is index i currently in the heap?
VSPINT MaxHeap_Contains_v2
(
    const max_heap_v2 *h,
    VSPINT             i
)
{
    return (h->lookup[i] != EMPTY) ;
}

/* On input:
   The child heaps of node p are 100% valid, but the value in node p may be
   out of order. lookup is valid for every node except possibly p.

   On output:
   The heap starting at node p is valid, and lookup is valid everywhere.
*/
static void MaxHeap_Bubble_Down_v2
(
    max_heap_v2 *h,
    VSPINT       p
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT    *heap = h->heap ;
    VSPINT  *lookup = h->lookup ;
    VSPINT       *c = h->c ;
    VSPINT numNodes = h->numNodes ;

    VSPINT Xe, Xleft, Xright ;
    VSPINT left, right, e, hleft, hright ;

    e  = heap[p] ;
    Xe = c[e] ;

    while ( TRUE )
    {
        left  = p * 2 ;
        right = left + 1 ;

        if ( right <= numNodes )
        {
            /* both left and right children of p are in the heap */
            hleft  = heap[left] ;
            hright = heap[right] ;
            Xleft  = c[hleft] ;
            Xright = c[hright] ;
            if (Xleft > Xright)
            {
                if (Xe < Xleft)
                {
                    heap[p]       = hleft ;
                    lookup[hleft] = p ;
                    p = left ;
                }
                else
                {
                    heap[p]   = e ;
                    lookup[e] = p ;
                    return ;
                }
            }
            else
            {
                if (Xe < Xright)
                {
                    heap[p]        = hright ;
                    lookup[hright] = p ;
                    p = right ;
                }
                else
                {
                    heap[p]   = e ;
                    lookup[e] = p ;
                    return ;
                }
            }
        }
        else
        {
            /* only the left child of p is in the heap */
            if ( left <= numNodes )
            {
                hleft = heap[left] ;
                Xleft = c[hleft] ;
                if (Xe < Xleft)
                {
                    heap[p]       = hleft ;
                    lookup[hleft] = p ;
                    p = left ;
                }
            }
            heap[p]   = e ;
            lookup[e] = p ;
            return ;
        }
    }
}

// Input  : h->heap[1..h->numNodes] is populated;
//          lookup[i] == EMPTY for all i NOT in the heap
// Output : heap satisfies the max-heap property; lookup valid everywhere
void MaxHeap_Heapify_v2
(
    max_heap_v2 *h
)
{
    for (VSPINT i = 1 ; i <= h->numNodes ; i++)
    {
        h->lookup[ h->heap[i] ] = i ;
    }

    for (VSPINT p = h->numNodes / 2 ; p >= 1 ; p--)
    {
        MaxHeap_Bubble_Down_v2(h, p) ;
    }

    #ifdef DEBUG
        MaxHeap_Check_v2(h) ;
    #endif
    return ;
}

// Bulk-load: discard current contents, load numElems indices, heapify. O(n).
void MaxHeap_Build_v2
(
    max_heap_v2  *h,
    const VSPINT *elems,
    VSPINT        numElems
)
{
    MaxHeap_Reset_v2(h) ;
    for (VSPINT i = 1 ; i <= numElems ; i++)
    {
        h->heap[i] = elems[i - 1] ;
    }
    h->numNodes = numElems ;
    MaxHeap_Heapify_v2(h) ;
    return ;
}

// ASSUMED: toInsert is NOT in the heap, and 0 <= toInsert < h->n
void MaxHeap_Insert_v2
(
    max_heap_v2 *h,
    VSPINT       toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        if (toInsert < 0 || toInsert >= h->n)
        {
            Error("Element out of range", DBG_INFO) ;
        }
        if (h->lookup[toInsert] != EMPTY)
        {
            Error("Element already in heap", DBG_INFO) ;
        }
        MaxHeap_Check_v2(h) ;
    #endif

    VSPINT    *heap = h->heap ;
    VSPINT  *lookup = h->lookup ;
    VSPINT       *c = h->c ;

    VSPINT node_par ;
    VSPINT node_child = ++(h->numNodes) ;
    VSPINT c_child    = c[toInsert] ;

    /* toInsert is perpetually considered the child, until its final place
       is found */
    while (node_child > 1)
    {
        node_par = node_child / 2 ;
        if ( c[ heap[node_par] ] < c_child )
        {
            heap[node_child]           = heap[node_par] ;
            lookup[ heap[node_child] ] = node_child ;
            node_child                 = node_par ;
        }
        else
        {
            break ;
        }
    }

    heap[node_child] = toInsert ;
    lookup[toInsert] = node_child ;

    #ifdef DEBUG
        MaxHeap_Check_v2(h) ;
    #endif
    return ;
}

// Remove the element at heap position `address`
void MaxHeap_Delete_v2
(
    max_heap_v2 *h,
    VSPINT       address
)
{
    #ifdef DEBUG
        MaxHeap_Check_v2(h) ;
    #endif

    if (address < 1 || address > h->numNodes)
    {
        return ;
    }

    VSPINT   *heap = h->heap ;
    VSPINT *lookup = h->lookup ;
    VSPINT      *c = h->c ;

    lookup[ heap[address] ] = EMPTY ; /* this element is leaving the heap */

    if (address == h->numNodes) /* element to delete is the last node */
    {
        h->numNodes-- ;
        #ifdef DEBUG
            MaxHeap_Check_v2(h) ;
        #endif
        return ;
    }

    VSPINT replacement = heap[ h->numNodes-- ] ;

    /* Determine if we need to bubble the replacement up or down */

    if ( address == 1 || c[ heap[address / 2] ] >= c[replacement] )
    {
        /* Bubble down (Bubble_Down_v2 sets lookup[replacement]) */
        heap[address] = replacement ;
        MaxHeap_Bubble_Down_v2(h, address) ;
    }
    else
    {
        /* Bubble up */
        while (address > 1)
        {
            VSPINT par = address / 2 ;
            if (c[ heap[par] ] < c[replacement])
            {
                heap[address]           = heap[par] ;
                lookup[ heap[address] ] = address ;
                address                 = par ;
            }
            else
            {
                break ;
            }
        }
        heap[address]       = replacement ;
        lookup[replacement] = address ;
    }

    #ifdef DEBUG
        MaxHeap_Check_v2(h) ;
    #endif
    return ;
}

// Delete an arbitrary element `toDel` (an index into c) from the heap.
// Does nothing if toDel is not in the heap.
void MaxHeap_Delete_Val_v2
(
    max_heap_v2 *h,
    VSPINT       toDel
)
{
    VSPINT address = h->lookup[toDel] ;
    if (address == EMPTY)
    {
        return ;
    }
    MaxHeap_Delete_v2(h, address) ;
    return ;
}

// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Val_v2
(
    const max_heap_v2 *h,
    VSPINT            *topVal
)
{
    *topVal = (h->numNodes != 0) ? h->heap[1] : EMPTY ;
    return ;
}

// Returns EMPTY if the heap is empty
void MaxHeap_Get_Top_Val_v2
(
    max_heap_v2 *h,
    VSPINT      *topVal
)
{
    if (h->numNodes == 0)
    {
        *topVal = EMPTY ;
        return ;
    }
    *topVal = h->heap[1] ;
    MaxHeap_Delete_v2(h, 1) ;
    return ;
}

// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Key_v2
(
    const max_heap_v2 *h,
    VSPINT            *topKey
)
{
    *topKey = (h->numNodes != 0) ? h->c[ h->heap[1] ] : EMPTY ;
    return ;
}

// Returns EMPTY if the heap is empty
void MaxHeap_Get_Top_Key_v2
(
    max_heap_v2 *h,
    VSPINT      *topKey
)
{
    if (h->numNodes == 0)
    {
        *topKey = EMPTY ;
        return ;
    }
    *topKey = h->c[ h->heap[1] ] ;
    MaxHeap_Delete_v2(h, 1) ;
    return ;
}

void MaxHeap_Check_v2
(
    const max_heap_v2 *h
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    /* Forward direction: every heap slot has a correct lookup entry */
    for (VSPINT i = 1 ; i <= h->numNodes ; i++)
    {
        VSPINT e = h->heap[i] ;
        if (e < 0 || e >= h->n)
        {
            ERR("Heap element out of range") ;
        }
        if (h->lookup[e] != i)
        {
            ERR("lookup inconsistent with heap") ;
        }
        if (i > 1)
        {
            VSPINT parent = i / 2 ;
            if (h->c[ h->heap[parent] ] < h->c[e])
            {
                ERR("Heap property violated") ;
            }
        }
    }

    /* Reverse direction: every non-EMPTY lookup entry points back */
    for (VSPINT i = 0 ; i < h->n ; i++)
    {
        VSPINT pos = h->lookup[i] ;
        if (pos == EMPTY) continue ;
        if (pos < 1 || pos > h->numNodes || h->heap[pos] != i)
        {
            ERR("Stale lookup entry") ;
        }
    }
}

/*---------------------------------------------------------------------------*/
// Min heap
/*---------------------------------------------------------------------------*/

min_heap_v2 *MinHeap_Init_v2
(
    VSPINT *c,
    VSPINT  n
)
{
    min_heap_v2 *h = malloc(sizeof(min_heap_v2)) ;
    if (h == NULL)
    {
        return NULL ;
    }

    h->heap   = malloc((size_t)(n + 1) * sizeof(VSPINT)) ;
    h->lookup = malloc((size_t) n      * sizeof(VSPINT)) ;
    if (h->heap == NULL || h->lookup == NULL)
    {
        free(h->heap) ;
        free(h->lookup) ;
        free(h) ;
        return NULL ;
    }

    h->c        = c ;
    h->n        = n ;
    h->numNodes = 0 ;

    for (VSPINT i = 0 ; i < n ; i++)
    {
        h->lookup[i] = EMPTY ;
    }
    return h ;
}

void MinHeap_Free_v2
(
    min_heap_v2 *h
)
{
    if (h == NULL)
    {
        return ;
    }
    free(h->heap) ;
    free(h->lookup) ;
    free(h) ;
    return ;
}

void MinHeap_Reset_v2
(
    min_heap_v2 *h
)
{
    h->numNodes = 0 ;
    for (VSPINT i = 0 ; i < h->n ; i++)
    {
        h->lookup[i] = EMPTY ;
    }
    return ;
}

VSPINT MinHeap_Contains_v2
(
    const min_heap_v2 *h,
    VSPINT             i
)
{
    return (h->lookup[i] != EMPTY) ;
}

static void MinHeap_Bubble_Down_v2
(
    min_heap_v2 *h,
    VSPINT       p
)
{
    VSPINT    *heap = h->heap ;
    VSPINT  *lookup = h->lookup ;
    VSPINT       *c = h->c ;
    VSPINT numNodes = h->numNodes ;

    VSPINT e, Xe, hleft, hright, Xleft, Xright ;
    VSPINT left, right ;

    e  = heap[p] ;
    Xe = c[e] ;

    while (TRUE)
    {
        left  = p << 1 ;     // p * 2
        right = left + 1 ;

        if (right <= numNodes)
        {
            hleft  = heap[left] ;
            hright = heap[right] ;
            Xleft  = c[hleft] ;
            Xright = c[hright] ;

            if (Xleft < Xright)
            {
                if (Xe > Xleft)
                {
                    heap[p]       = hleft ;
                    lookup[hleft] = p ;
                    p = left ;
                }
                else
                {
                    heap[p]   = e ;
                    lookup[e] = p ;
                    return ;
                }
            }
            else
            {
                if (Xe > Xright)
                {
                    heap[p]        = hright ;
                    lookup[hright] = p ;
                    p = right ;
                }
                else
                {
                    heap[p]   = e ;
                    lookup[e] = p ;
                    return ;
                }
            }
        }
        else
        {
            if (left <= numNodes)
            {
                hleft = heap[left] ;
                Xleft = c[hleft] ;
                if (Xe > Xleft)
                {
                    heap[p]       = hleft ;
                    lookup[hleft] = p ;
                    p = left ;
                }
            }
            heap[p]   = e ;
            lookup[e] = p ;
            return ;
        }
    }
}

void MinHeap_Heapify_v2
(
    min_heap_v2 *h
)
{
    for (VSPINT i = 1 ; i <= h->numNodes ; i++)
    {
        h->lookup[ h->heap[i] ] = i ;
    }

    for (VSPINT p = h->numNodes / 2 ; p >= 1 ; p--)
    {
        MinHeap_Bubble_Down_v2(h, p) ;
    }
    #ifdef DEBUG
        MinHeap_Check_v2(h) ;
    #endif
}

void MinHeap_Build_v2
(
    min_heap_v2  *h,
    const VSPINT *elems,
    VSPINT        numElems
)
{
    MinHeap_Reset_v2(h) ;
    for (VSPINT i = 1 ; i <= numElems ; i++)
    {
        h->heap[i] = elems[i - 1] ;
    }
    h->numNodes = numElems ;
    MinHeap_Heapify_v2(h) ;
    return ;
}

void MinHeap_Insert_v2
(
    min_heap_v2 *h,
    VSPINT       toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        MinHeap_Check_v2(h) ;
        if (toInsert < 0 || toInsert >= h->n)
        {
            Error("Element out of range", DBG_INFO) ;
        }
        if (h->lookup[toInsert] != EMPTY)
        {
            Error("Element already in heap", DBG_INFO) ;
        }
    #endif

    VSPINT    *heap = h->heap ;
    VSPINT  *lookup = h->lookup ;
    VSPINT       *c = h->c ;

    VSPINT node_child = ++(h->numNodes) ;
    VSPINT c_child    = c[toInsert] ;
    VSPINT node_par ;

    while (node_child > 1)
    {
        node_par = node_child >> 1 ; // node_child / 2
        if (c[heap[node_par]] > c_child)
        {
            heap[node_child]         = heap[node_par] ;
            lookup[heap[node_child]] = node_child ;
            node_child               = node_par ;
        }
        else
        {
            break ;
        }
    }

    heap[node_child] = toInsert ;
    lookup[toInsert] = node_child ;

    #ifdef DEBUG
        MinHeap_Check_v2(h) ;
    #endif
}

// Delete the element at heap position `address`
void MinHeap_Delete_v2
(
    min_heap_v2 *h,
    VSPINT       address
)
{
    #ifdef DEBUG
        MinHeap_Check_v2(h) ;
    #endif

    if (address < 1 || address > h->numNodes) return ;

    VSPINT   *heap = h->heap ;
    VSPINT *lookup = h->lookup ;
    VSPINT      *c = h->c ;

    lookup[heap[address]] = EMPTY ; /* this element is leaving the heap */

    if (address == h->numNodes)
    {
        h->numNodes-- ;
        #ifdef DEBUG
            MinHeap_Check_v2(h) ;
        #endif
        return ;
    }

    VSPINT replacement = heap[ h->numNodes-- ] ;

    if (address == 1 || c[heap[address >> 1]] <= c[replacement])
    {
        /* Bubble down (Bubble_Down_v2 sets lookup[replacement]) */
        heap[address] = replacement ;
        MinHeap_Bubble_Down_v2(h, address) ;
    }
    else
    {
        while (address > 1)
        {
            VSPINT par = address >> 1 ;
            if (c[heap[par]] > c[replacement])
            {
                heap[address]         = heap[par] ;
                lookup[heap[address]] = address ;
                address               = par ;
            }
            else
            {
                break ;
            }
        }
        heap[address]       = replacement ;
        lookup[replacement] = address ;
    }
    #ifdef DEBUG
        MinHeap_Check_v2(h) ;
    #endif
}

// Delete an arbitrary element `toDel` (an index into c) from the heap.
// Does nothing if toDel is not in the heap.
void MinHeap_Delete_Val_v2
(
    min_heap_v2 *h,
    VSPINT       toDel
)
{
    VSPINT address = h->lookup[toDel] ;
    if (address == EMPTY) return ;
    MinHeap_Delete_v2(h, address) ;
}

void MinHeap_See_Top_Val_v2
(
    const min_heap_v2 *h,
    VSPINT            *topVal
)
{
    *topVal = (h->numNodes != 0) ? h->heap[1] : EMPTY ;
}

void MinHeap_Get_Top_Val_v2
(
    min_heap_v2 *h,
    VSPINT      *topVal
)
{
    if (h->numNodes == 0)
    {
        *topVal = EMPTY ;
        return ;
    }
    *topVal = h->heap[1] ;
    MinHeap_Delete_v2(h, 1) ;
}

void MinHeap_See_Top_Key_v2
(
    const min_heap_v2 *h,
    VSPINT            *topKey
)
{
    *topKey = (h->numNodes != 0) ? h->c[ h->heap[1] ] : EMPTY ;
}

void MinHeap_Get_Top_Key_v2
(
    min_heap_v2 *h,
    VSPINT      *topKey
)
{
    if (h->numNodes == 0)
    {
        *topKey = EMPTY ;
        return ;
    }
    *topKey = h->c[ h->heap[1] ] ;
    MinHeap_Delete_v2(h, 1) ;
}

void MinHeap_Check_v2
(
    const min_heap_v2 *h
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    for (VSPINT i = 1 ; i <= h->numNodes ; i++)
    {
        VSPINT e = h->heap[i] ;
        if (e < 0 || e >= h->n)
        {
            ERR("Heap element out of range") ;
        }
        if (h->lookup[e] != i)
        {
            ERR("lookup inconsistent with heap") ;
        }
        if (i > 1)
        {
            VSPINT parent = i / 2 ;
            if (h->c[ h->heap[parent] ] > h->c[e])
            {
                ERR("Heap property violated") ;
            }
        }
    }

    for (VSPINT i = 0 ; i < h->n ; i++)
    {
        VSPINT pos = h->lookup[i] ;
        if (pos == EMPTY) continue ;
        if (pos < 1 || pos > h->numNodes || h->heap[pos] != i)
        {
            ERR("Stale lookup entry") ;
        }
    }
}
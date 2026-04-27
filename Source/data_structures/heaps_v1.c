/******************************************************************************
VERSION 0 (lookup-free)
Simple min/max heap functions. The heap stores indices into the c array;
heap keys are the c values at those indices.

API CHANGE FROM ORIGINAL:
    *_Delete_v1 now takes a heap position (address) rather than an element
    value, because without the lookup array there is no O(1) way to translate
    a value to its heap position. The Get_Top_* helpers always pass address=1.
*******************************************************************************/

#include "core.h"
#include "data_structures.h"
#include "utils.h"

/*----------------------------------------------------------------------------*/
// Max heap
/*----------------------------------------------------------------------------*/

/* On input:
   The child heaps of node p are 100% valid, but the value in node p may be
   out of order.

   On output:
   The heap starting at node p is valid
*/
static void MaxHeap_Bubble_Down_v1
(
    VSPINT       *heap,
    VSPINT    numNodes, // Number of nodes currently in the heap
    VSPINT          *c,
    VSPINT           p  // start at node p in the heap
)
{
    SET_FUNCTION_COUNTER ;

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
                /* traverse down the left heap */
                if (Xe < Xleft)
                {
                    heap[p] = hleft ;
                    p = left ;
                }
                else
                {
                    heap[p] = e ;
                    return ;
                }
            }
            else
            {
                /* traverse down the right heap */
                if (Xe < Xright)
                {
                    heap[p] = hright ;
                    p = right ;
                }
                else
                {
                    heap[p] = e ;
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
                    heap[p] = hleft ;
                    p = left ;
                }
            }
            heap[p] = e ;
            return ;
        }
    }
}

// Input  : heap array is populated
// Output : heap array satisfies the max-heap property
void MaxHeap_Heapify_v1
(
    VSPINT       *heap,
    VSPINT    numNodes,
    VSPINT          *c
)
{
    for (VSPINT p = numNodes / 2 ; p >= 1 ; p--)
    {
        MaxHeap_Bubble_Down_v1(heap, numNodes, c, p) ;
    }

    #ifdef DEBUG
        MaxHeap_Check_v1(heap, numNodes, c) ;
    #endif
    return ;
}

// ASSUMED: toInsert is NOT in the queue
void MaxHeap_Insert_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT  toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        if(Is_Array_Member(toInsert, heap, 1, *numNodes + 1))
        {
            Error("Element already in heap", DBG_INFO) ;
        }
        MaxHeap_Check_v1(heap, *numNodes, c) ;
    #endif

    VSPINT node_par ;
    VSPINT node_child = ++(*numNodes) ;
    VSPINT c_child    = c[toInsert] ;

    /* toInsert is perpetually considered the child, until its final place
       is found */
    while(node_child > 1)
    {
        node_par = node_child / 2 ;
        if( c[ heap[node_par] ] < c_child )
        {
            heap[node_child] = heap[node_par] ;
            node_child       = node_par ;
        }
        else
        {
            break ;
        }
    }

    heap[node_child] = toInsert ;

    #ifdef DEBUG
        MaxHeap_Check_v1(heap, *numNodes, c) ;
    #endif
    return ;
}

// Remove the element at heap position `address`
// (formerly took an element value `toDel`; without lookup we use position)
void MaxHeap_Delete_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   address
)
{
    #ifdef DEBUG
        MaxHeap_Check_v1(heap, *numNodes, c) ;
    #endif

    if (address < 1 || address > *numNodes)
    {
        return ;
    }
    if (address == *numNodes) /* element to delete is the last node */
    {
        (*numNodes)-- ;
        #ifdef DEBUG
            MaxHeap_Check_v1(heap, *numNodes, c) ;
        #endif
        return ;
    }

    VSPINT replacement = heap[ (*numNodes)-- ] ;

    /* Determine if we need to bubble the replacement up or down */

    if( address == 1 || c[ heap[address / 2] ] >= c[replacement] )
    {
        /* Bubble down */
        heap[address] = replacement ;
        MaxHeap_Bubble_Down_v1(heap, *numNodes, c, address) ;
    }
    else
    {
        /* Bubble up */
        while(address > 1)
        {
            VSPINT par = address / 2 ;
            if(c[ heap[par] ] < c[replacement])
            {
                heap[address] = heap[par] ;
                address = par ;
            }
            else
            {
                break ;
            }
        }
        heap[address] = replacement ;
    }

    #ifdef DEBUG
        MaxHeap_Check_v1(heap, *numNodes, c) ;
    #endif
    return ;
}

// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topVal
)
{
    *topVal = (numNodes != 0) ? heap[1] : EMPTY ;
    return ;
}

// Returns EMPTY if the heap is empty
void MaxHeap_Get_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topVal
)
{
    if( *numNodes == 0 )
    {
        *topVal = EMPTY ;
        return ;
    }
    *topVal = heap[1] ;
    MaxHeap_Delete_v1(heap, numNodes, c, 1) ;
    return ;
}

// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topKey
)
{
    *topKey = (numNodes != 0) ? c[heap[1]] : EMPTY ;
    return ;
}

// Returns EMPTY if the heap is empty
void MaxHeap_Get_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topKey
)
{
    if( *numNodes == 0 )
    {
        *topKey = EMPTY ;
        return ;
    }
    *topKey = c[ heap[1] ] ;
    MaxHeap_Delete_v1(heap, numNodes, c, 1) ;
    return ;
}

void MaxHeap_Check_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    for (VSPINT i = 1 ; i <= numNodes ; i++)
    {
        VSPINT e = heap[i] ;
        if( e < 0 )
        {
            ERR("Negative value in heap") ;
        }
        // Check heap property
        if (i > 1)
        {
            VSPINT parent = i / 2 ;
            if (c[ heap[parent] ] < c[e])
            {
                ERR("Heap property violated") ;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// Min heap
/*---------------------------------------------------------------------------*/

void MinHeap_Bubble_Down_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c,
    VSPINT         p
)
{
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
                    heap[p] = hleft ;
                    p = left ;
                }
                else
                {
                    heap[p] = e ;
                    return ;
                }
            }
            else
            {
                if (Xe > Xright)
                {
                    heap[p] = hright ;
                    p = right ;
                }
                else
                {
                    heap[p] = e ;
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
                    heap[p] = hleft ;
                    p = left ;
                }
            }
            heap[p] = e ;
            return ;
        }
    }
}

void MinHeap_Heapify_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
)
{
    for (VSPINT p = numNodes / 2 ; p >= 1 ; p--)
    {
        MinHeap_Bubble_Down_v1(heap, numNodes, c, p) ;
    }
    #ifdef DEBUG
        MinHeap_Check_v1(heap, numNodes, c) ;
    #endif
}

void MinHeap_Insert_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT  toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        MinHeap_Check_v1(heap, *numNodes, c) ;
        if(Is_Array_Member(toInsert, heap, 1, *numNodes + 1))
        {
            Error("Element already in heap", DBG_INFO) ;
        }
    #endif

    VSPINT node_child = ++(*numNodes) ;
    VSPINT c_child    = c[toInsert] ;
    VSPINT node_par ;

    while (node_child > 1)
    {
        node_par = node_child >> 1 ; // node_child / 2
        if (c[heap[node_par]] > c_child)
        {
            heap[node_child] = heap[node_par] ;
            node_child       = node_par ;
        }
        else
        {
            break ;
        }
    }

    heap[node_child] = toInsert ;

    #ifdef DEBUG
        MinHeap_Check_v1(heap, *numNodes, c) ;
    #endif
}

// Delete the element at heap position `address`
void MinHeap_Delete_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   address
)
{
    #ifdef DEBUG
        MinHeap_Check_v1(heap, *numNodes, c) ;
    #endif

    if (address < 1 || address > *numNodes) return ;
    if (address == *numNodes)
    {
        (*numNodes)-- ;
        #ifdef DEBUG
            MinHeap_Check_v1(heap, *numNodes, c) ;
        #endif
        return ;
    }

    VSPINT replacement = heap[(*numNodes)--] ;

    if (address == 1 || c[heap[address >> 1]] <= c[replacement])
    {
        heap[address] = replacement ;
        MinHeap_Bubble_Down_v1(heap, *numNodes, c, address) ;
    }
    else
    {
        while (address > 1)
        {
            VSPINT par = address >> 1 ;
            if (c[heap[par]] > c[replacement])
            {
                heap[address] = heap[par] ;
                address       = par ;
            }
            else
            {
                break ;
            }
        }
        heap[address] = replacement ;
    }
    #ifdef DEBUG
        MinHeap_Check_v1(heap, *numNodes, c) ;
    #endif
}

void MinHeap_See_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topVal
)
{
    *topVal = (numNodes != 0) ? heap[1] : EMPTY ;
}

void MinHeap_Get_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topVal
)
{
    if (*numNodes == 0)
    {
        *topVal = EMPTY ;
        return ;
    }
    *topVal = heap[1] ;
    MinHeap_Delete_v1(heap, numNodes, c, 1) ;
}

void MinHeap_See_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topKey
)
{
    *topKey = (numNodes != 0) ? c[heap[1]] : EMPTY ;
}

void MinHeap_Get_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topKey
)
{
    if (*numNodes == 0)
    {
        *topKey = EMPTY ;
        return ;
    }
    *topKey = c[heap[1]] ;
    MinHeap_Delete_v1(heap, numNodes, c, 1) ;
}

void MinHeap_Check_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    for (VSPINT i = 1 ; i <= numNodes ; i++)
    {
        VSPINT e = heap[i] ;
        if (e < 0)
        {
            ERR("Negative value in heap") ;
        }
        if (i > 1)
        {
            VSPINT parent = i / 2 ;
            if (c[ heap[parent] ] > c[e])
            {
                ERR("Heap property violated") ;
            }
        }
    }
}
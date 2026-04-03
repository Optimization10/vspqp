#include "core.h"
/******************************************************************************
VERSION 0
These are simple min/max heap functions that are designed to be flexible by
taking in arrays, and not any particular data struct. The heap values are
indices into the c array, and the heap keys are the c_i values

ASSUMED: The pointer array may be shared between multiple heaps. Hence,
if i is in the heap, then
    ptrs[i] = node number where i is located
However, if i is not in the heap, there is no guarantee on ptrs[i]
*******************************************************************************/

#include "shared_api.h"
#include "ds_api.h"


/*----------------------------------------------------------------------------*/
// Max heap
/*----------------------------------------------------------------------------*/

// Input  : heap array is populated and lookup is correct
// Output : heap array is a heap and lookup is correct
void MaxHeap_Heapify_v0
(
    VSPINT           *heap,
    VSPINT         *lookup,
    VSPINT        numNodes, // Number of nodes currently in the heap
    VSPINT              *c
)
{
    for (VSPINT p = numNodes / 2  ; p >= 1  ; p--)
    {
        MaxHeap_Bubble_Down_v0(heap, lookup, numNodes, c, p) ;
    }

    #ifdef VSP_DEBUG
        MaxHeap_Check_v0(heap, lookup, numNodes, c) ;
    #endif
    return ;
}

/* On input:
   The child heaps of node p are 100% valid, but the value in node p may be
   out of order.

   On output: 
   The heap starting at node p is valid
*/
void MaxHeap_Bubble_Down_v0
(
    VSPINT           *heap,
    VSPINT         *lookup,
    VSPINT        numNodes, // Number of nodes currently in the heap
    VSPINT              *c,
    VSPINT               p  // start at node p in the heap
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT Xe, Xleft, Xright ;
    VSPINT left, right, e, hleft, hright ;

    e = heap [p] ;
    Xe = c [e] ;

    while ( TRUE )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= numNodes )
        {
            /* both left and right children of p are in the heap */
            hleft = heap [left] ;
            hright = heap [right] ;
            Xleft = c [hleft] ;
            Xright = c [hright] ;
            if ((Xleft > Xright))
            {
                /* traverse down the left heap */
                if ((Xe < Xleft))
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    lookup [hleft] = p ;
                    p = left ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    lookup [e] = p ;
                    return ;
                }
            }
            else
            {
                /* traverse down the right heap */
                if ((Xe < Xright)) /* for min heap: (Xe > Xright) */
                {
                    /* swap e with the right child */
                    heap [p] = hright ;
                    lookup [hright] = p ;
                    p = right ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    lookup [e] = p ;
                    return ;
                }
            }
        }
        else
        {
            /* only the left child of p is in the heap */
            if ( left <= numNodes )
            {
                /* traverse down the left heap */
                hleft = heap [left] ;
                Xleft = c [hleft] ;
                if ((Xe < Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    lookup [hleft] = p ;
                    p = left ;
                }
            }
            /* element e, with value Xe, has found its final home */
            heap [p] = e ;
            lookup [e] = p ;
            return ;
        }
    }
}


// ASSUMED: to_insert is NOT in the queue
void MaxHeap_Insert_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT  *numNodes,
    VSPINT         *c,
    VSPINT   toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        if(Is_Array_Member(toInsert, heap, 1, *numNodes + 1) )
        {
            Error("Element already in heap", DBG_INFO) ;
        }

        MaxHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif

    VSPINT node_par ;
    VSPINT node_child = ++(*numNodes)  ;
    VSPINT c_child  = c[toInsert] ;

    /* to_insert is perpetually considered the child, until its final place
       is found */
    while(node_child > 1)
    {
        node_par = node_child / 2 ;
        if( c[ heap[node_par] ] < c_child )
        {
            // Write the parent to the child
            heap[node_child]         = heap[node_par] ;
            lookup[ heap[node_par] ] = node_child ;
            node_child = node_par ;
        }
        else
        {
            break ;
        }
    }

    
    heap[node_child] = toInsert ;
    lookup[toInsert] = node_child ;

    #ifdef VSP_DEBUG
        MaxHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif

    return ;
}



// Remove the value "toDel" from the queue
// ASSUMED: toDel is in the heap
void MaxHeap_Delete_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT  *numNodes,
    VSPINT         *c,
    VSPINT      toDel
)
{
    #ifdef VSP_DEBUG
        MaxHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif

    VSPINT address  = lookup[toDel] ;
    if(address == EMPTY)
    {
        return ;
    }
    else if(address == *numNodes) /* toDel is in the last node */
    {
        (*numNodes)-- ;
        #ifdef VSP_DEBUG
            MaxHeap_Check_v0(heap, lookup, *numNodes, c) ;
        #endif
        return ;
    }

    VSPINT replacement = heap[ (*numNodes)-- ] ;
    lookup[toDel] = EMPTY ;

    /* Determine if we need to bubble the replacement up or down */

    if( address == 1 || c[ heap[address / 2] ] >= c[replacement] )
    {
        /* Bubble down */
        heap[address] = replacement ;
        lookup[replacement] = address ;
        MaxHeap_Bubble_Down_v0(heap, lookup, *numNodes, c, address) ;
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
                lookup[heap[par]] = address ;
                address = par ;
            }
            else
            {
                break ;
            }
        }
        heap[address] = replacement ;
        lookup[replacement] = address ;
    }

    #ifdef VSP_DEBUG
        MaxHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif
    return ;
}


// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Val_v0
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
void MaxHeap_Get_Top_Val_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT  *numNodes,
    VSPINT         *c,
    VSPINT    *topVal
)
{
    if( *numNodes == 0 )
    {
        *topVal = EMPTY ;
        return ;
    }
    *topVal = heap[1] ;
    MaxHeap_Delete_v0(heap, lookup, numNodes, c, *topVal) ;
    return ;
}

// Returns EMPTY if heap is empty
void MaxHeap_See_Top_Key_v0
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
void MaxHeap_Get_Top_Key_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT  *numNodes,
    VSPINT         *c,
    VSPINT    *topKey
)
{
    if( *numNodes == 0 )
    {
        *topKey = EMPTY ;
        return ;
    }
    *topKey = c[ heap[1] ] ;
    MaxHeap_Delete_v0(heap, lookup, numNodes, c, heap[1]) ;
    return ;
}


void MaxHeap_Check_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT   numNodes,
    VSPINT         *c
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    for (VSPINT i = 1; i <= numNodes; i++)
    {
        VSPINT e = heap[i];
        if( e < 0 )
        {
            ERR("Negative value in heap") ;
        }

        // Check if lookup and heap are consistent
        if (lookup[e] != i)
        {
            ERR("Ptrs into heap incorrect") ;
        }
        // Check heap property
        if (i > 1)
        {
            VSPINT parent = i / 2;
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

void MinHeap_Bubble_Down_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT   numNodes,
    VSPINT  *c,
    VSPINT   p
)
{
    VSPINT e, Xe, hleft, hright, Xleft, Xright;
    VSPINT left, right;

    e  = heap[p];
    Xe = c[e];

    while (TRUE)
    {
        left  = p << 1;      // p * 2
        right = left + 1;

        if (right <= numNodes)
        {
            hleft  = heap[left];
            hright = heap[right];
            Xleft  = c[hleft];
            Xright = c[hright];

            if (Xleft < Xright)
            {
                if (Xe > Xleft)
                {
                    heap[p]     = hleft;
                    lookup[hleft] = p;
                    p  = left;
                }
                else
                {
                    heap[p]   = e;
                    lookup[e] = p;
                    return;
                }
            }
            else
            {
                if (Xe > Xright)
                {
                    heap[p]      = hright;
                    lookup[hright] = p;
                    p   = right;
                }
                else
                {
                    heap[p]   = e;
                    lookup[e] = p;
                    return;
                }
            }
        }
        else
        {
            if (left <= numNodes)
            {
                hleft = heap[left];
                Xleft = c[hleft];
                if (Xe > Xleft)
                {
                    heap[p]     = hleft;
                    lookup[hleft] = p;
                    p   = left;
                }
            }
            heap[p]   = e;
            lookup[e] = p;
            return;
        }
    }
}

void MinHeap_Heapify_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT   numNodes,
    VSPINT  *c
)
{
    for (VSPINT p = numNodes / 2; p >= 1; p--)
    {
        MinHeap_Bubble_Down_v0(heap, lookup, numNodes, c, p);
    }
    #ifdef VSP_DEBUG
        MinHeap_Check_v0(heap, lookup, numNodes, c) ;
    #endif
}

void MinHeap_Insert_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toInsert
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        MinHeap_Check_v0(heap, lookup, *numNodes, c) ;

        if(Is_Array_Member(toInsert, heap, 1, *numNodes + 1) )
        {
            Error("Element already in heap", DBG_INFO) ;
        }
    #endif

    VSPINT node_child = ++(*numNodes);
    VSPINT c_child    = c[toInsert];
    VSPINT node_par;

    while (node_child > 1)
    {
        node_par = node_child >> 1; // node_child / 2
        if (c[heap[node_par]] > c_child)
        {
            heap[node_child]           = heap[node_par];
            lookup[heap[node_par]]     = node_child;
            node_child                 = node_par;
        }
        else
        {
            break;
        }
    }

    heap[node_child]   = toInsert;
    lookup[toInsert]   = node_child;

    #ifdef VSP_DEBUG
        MinHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif
}

// ASSUMED: toDel is in the heap
void MinHeap_Delete_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toDel
)
{
    #ifdef VSP_DEBUG
        MinHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif

    VSPINT address = lookup[toDel];
    if (address == EMPTY) return;
    if (address == *numNodes)
    {
        (*numNodes)--;
        #ifdef VSP_DEBUG
            MinHeap_Check_v0(heap, lookup, *numNodes, c) ;
        #endif
        return;
    }

    VSPINT replacement = heap[(*numNodes)--];
    lookup[toDel]      = EMPTY;

    if (address == 1 || c[heap[address >> 1]] <= c[replacement])
    {
        heap[address]    = replacement;
        lookup[replacement] = address;
        MinHeap_Bubble_Down_v0(heap, lookup, *numNodes, c, address);
    }
    else
    {
        while (address > 1)
        {
            VSPINT par = address >> 1;
            if (c[heap[par]] > c[replacement])
            {
                heap[address]           = heap[par];
                lookup[heap[par]]       = address;
                address                 = par;
            }
            else
            {
                break;
            }
        }
        heap[address]     = replacement;
        lookup[replacement] = address;
    }
    #ifdef VSP_DEBUG
        MinHeap_Check_v0(heap, lookup, *numNodes, c) ;
    #endif
}

void MinHeap_See_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topVal
)
{
    *topVal = (numNodes != 0) ? heap[1] : EMPTY;
}

void MinHeap_Get_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topVal
)
{
    if (*numNodes == 0)
    {
        *topVal = EMPTY;
        return;
    }
    *topVal = heap[1];
    MinHeap_Delete_v0(heap, lookup, numNodes, c, *topVal);
}

void MinHeap_See_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topKey
)
{
    *topKey = (numNodes != 0) ? c[heap[1]] : EMPTY;
}

void MinHeap_Get_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topKey
)
{
    if (*numNodes == 0)
    {
        *topKey = EMPTY;
        return;
    }
    *topKey = c[heap[1]];
    MinHeap_Delete_v0(heap, lookup, numNodes, c, heap[1]);
}


void MinHeap_Check_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT   numNodes,
    VSPINT  *c
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE;

    for (VSPINT i = 1; i <= numNodes; i++)
    {
        VSPINT e = heap[i];
        if (e < 0)
        {
            ERR("Negative value in heap");
        }

        // Check if lookup and heap are consistent
        if (lookup[e] != i)
        {
            ERR("Ptrs into heap incorrect");
        }

        // Check min-heap property
        if (i > 1)
        {
            VSPINT parent = i / 2;
            if (c[ heap[parent] ] > c[e])
            {
                ERR("Heap property violated");
            }
        }
    }
}
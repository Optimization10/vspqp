#include "ds_api.h"

#define INSERT_IN_LIST(newListMem, insertLoc, list, pointers) \
                        (list)[ (insertLoc) ] = (newListMem) ; \
                        (pointers)[ (newListMem) ] = (insertLoc)


// For these heaps, the values and keys are the same thing
// We are NOT heapifing indices of some array by the array values
// We ARE heapifying values directly


/*----------------------------------------------------------------------------*/
/*   MIN HEAP FUNCTIONS                                                       */
/*----------------------------------------------------------------------------*/
// Constructor for Min_Heap_v4
Min_Heap_v4 *Min_Heap_Init_v4
(
    VSPINT gamma,
    VSPINT    *c,
    VSPINT     n  // Number of vertices at current level
)
{
    Min_Heap_v4* Heap = malloc(sizeof(*Heap)) ;
    Heap->c = c ;

    Heap->heap        = malloc( ( gamma + 1 ) * sizeof(VSPINT)) ;
    Heap->heapPtrs    = malloc( ( gamma + 1 ) * sizeof(VSPINT)) ;
    for(VSPINT i = 0 ; i < gamma + 1 ; i++)
    {
        Heap->heapPtrs[i] = EMPTY ;
    }
    Heap->sizeHeap = 0 ;
    Heap->gamma = gamma ;
    Heap->n     = n ;

    #ifdef VSP_DEBUG
        Min_Heap_Check_v4(Heap) ;
    #endif

    return Heap;
}

// Destructor for Min_Heap_v4
void Min_Heap_Destroy_v4(Min_Heap_v4 *heap)
{
    free(heap->heap);
    free(heap->heapPtrs);
    free(heap);
}

void Min_Heap_Check_v4
(
    Min_Heap_v4 *Heap
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT n = Heap->n ;
    VSPINT i, p, cval ;
    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;
    VSPINT  sizeHeap = Heap->sizeHeap ;
    VSPINT  gamma = Heap->gamma ;
    VSPINT *c = Heap->c ;
    VSP_ASSERT(n > 0) ;

    // Make sure size does not exceed maximum
    if (sizeHeap > gamma)
    {
        ERR("Heap size exceeds gamma") ;
    }

    /*------------------------------------------------------------------------*/
    // Check 
    // 1) Each value in the heap is in the correct range and it a c value
    // 2) The heap property is satisfies
    // 3) The pointer array points to it
    /*------------------------------------------------------------------------*/
    for(VSPINT par = 1 ; par <= sizeHeap ; par++)
    {
        VSPINT left   = 2 * par ;
        VSPINT right  = 2 * par + 1 ;
        VSPINT cpar   = heap[par] ;

        // Check that cpar is a value in the c array
        if( !Is_Array_Member(cpar, c, 0, n) )
        {
            ERR("Invalid entry into array") ;
        }
        
        // Check that the entry is in the valid range
        if( !IN_RANGE(cpar, 1, gamma + 1) )
        {
            ERR("Invalid entry in heap") ;
        }

        // Check the heap property
        if(left <= sizeHeap)
        {
            VSPINT cleft  = heap[left] ;
            if( !( cleft >= cpar ) )
            {
                ERR("Heap Property Violated") ;
            }
        }

        if(right <= sizeHeap)
        {
            VSPINT cright = heap[right] ;
            if( !( cright >= cpar ) )
            {
                ERR("Heap Property Violated") ;
            }
        }

        // Check the pointer
        if( ptrs[cpar] != par )
        {
            ERR("Heap pointers wrong") ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Check
    // 1) ptrs[cval] == -1 if and only if cval is not in the heap 
    /*------------------------------------------------------------------------*/
    for(VSPINT cval = 1 ; cval <= gamma ; cval++)
    {
        if( ptrs[cval] == EMPTY )
        {
            if( Is_Array_Member(cval, heap, 1, sizeHeap + 1) )
            {
                ERR("cval in heap but pointer says it's not") ;
            }
        }
        else
        {
            if( !Is_Array_Member(cval, heap, 1, sizeHeap + 1) )
            {
                ERR("Pointers say cval is in heap but it's actually not") ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Check the heap has no duplicated
    /*------------------------------------------------------------------------*/
    VSPINT *cvalsAppearing = calloc((gamma + 1), sizeof(VSPINT) ) ;
    for(p = 1 ; p <= sizeHeap ; p++)
    {
        cval = heap[p] ;
        cvalsAppearing[cval]++ ;
    }

    for(i = 0 ; i <= gamma ; i++)
    {
        if(cvalsAppearing[i] != 0 && cvalsAppearing[i] != 1)
        {
            ERR("Duplicate in heap") ;
        }
    }
    free(cvalsAppearing) ;
}

/* On input:
The subtree of "heap" starting at node with index "p" is a heap,
except that it's root may be out of order.

On output: 
subtree of "heap" starting at node with index "p" will be a subheap. */
void Min_Heap_Bubble_Down_v4
(
    Min_Heap_v4     *Heap,
    VSPINT              p  // start at node p in the heap
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT left, right, e, hleft, hright ;

    VSPINT *heap  = Heap->heap ;
    VSPINT *ptrs  = Heap->heapPtrs ;
    VSPINT  nheap = Heap->sizeHeap ;

    e  = heap[p] ;
    while ( 1 )
    {
        left  = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            // Both left and right children of p are in the heap
            hleft  = heap [left] ;
            hright = heap [right] ;

            // We MUST swap with the smaller child, there is no choice. If we
            // swap with the larger child, then the new "root" will violate
            // the heap property because it will be larger than one of its
            // children
            if ( hleft < hright )
            {
                // Traverse down the left heap
                if ( e > hleft )
                {
                    // Swap e with the left child
                    INSERT_IN_LIST(hleft, p, heap, ptrs) ;
                    p = left ;
                }
                else
                {
                    // Element e, with value ce, has found its final home
                    INSERT_IN_LIST(e, p, heap, ptrs) ;
                    return ;
                }
            }
            else
            {
                // Traverse down the right heap
                if ( e > hright )
                {
                    INSERT_IN_LIST(hright, p, heap, ptrs) ;
                    p = right ;
                }
                else
                {
                    // Element e, with value ce, has found its final home
                    INSERT_IN_LIST(e, p, heap, ptrs) ;
                    return ;
                }
            }
        }
        else
        {
            // Only the left child of p is in the heap
            if ( left <= nheap )
            {
                // Traverse down the left heap
                hleft = heap [left] ;
                if ( e > hleft )
                {
                    // Swap e with the left child
                    INSERT_IN_LIST(hleft, p, heap, ptrs) ;
                    p = left ;
                }
            }
            // Element e, with value ce, has found its final home
            INSERT_IN_LIST(e, p, heap, ptrs) ;
            return ;
        }
    }
}


// INPUT: The past lineage of p is all in order, i.e. 
// parent >= grandparent >= grand-grandparent .... 
// but p MAY be < parent
//
// Any children of p is ignored. They may or may not satisfy the heap property
//
// OUTPUT: the value residing at node p is in the correct place in its lineage
void Min_Heap_Bubble_Up_v4
(
    Min_Heap_v4         *Heap,
    VSPINT                  p, // node in the heap to bubble up
    VSPINT        valToBubble // value (equivalently, key) stored initially in node p
)
{
    VSPINT nodePar, valPar, keyPar ;
    VSPINT *heap  = Heap->heap ;
    VSPINT *ptrs  = Heap->heapPtrs ;
    VSPINT  n     = Heap->sizeHeap ;

    
    while(p > 1)
    {
        nodePar = p / 2 ;
        valPar  = heap[nodePar] ;
        if( valPar > valToBubble )
        {
            // Overwrite the child node with the parent node
            INSERT_IN_LIST(valPar, p, heap, ptrs) ;
            p = nodePar ; // TODO: Consider p = p / 2
        }
        else
        {
            break ;
        }
    }

    // Write valToBubble into to its proper location in the heap
    INSERT_IN_LIST(valToBubble, p, heap, ptrs) ;
    return ;
}


// Takes the array "heap" and heapifies it
void Min_Heap_Build_v4
(
    Min_Heap_v4 *Heap
)
{
    VSPINT n = Heap->sizeHeap ;
    for (VSPINT p = n / 2  ; p >= 1  ; p--)
    {
        Min_Heap_Bubble_Down_v4(Heap, p)  ;
    }
}


// Note: This function DOES make sure that a duplicate is not inserted.
// So before inserting a value, it checks that it is not already in the heap
// Trick to make this more efficient: when we bubble up and swap two nodes,
// it is not neccesary to actually write the child node to the parent node.
// It is sufficient to write the parent node to the child node. At the very end,
// when i can no longer bubble up, it can be written to the current child
void Min_Heap_Insert_v4
(
    Min_Heap_v4         *Heap,
    VSPINT        valToInsert  // value/hey to insert
)
{
    #ifdef VSP_DEBUG
        Min_Heap_Check_v4(Heap) ;
    #endif
    VSPINT *heap  = Heap->heap ;
    VSPINT *ptrs  = Heap->heapPtrs ;

    if(ptrs[valToInsert] != EMPTY)
    {
        return ;
    }

    // Insert the value at the end and bubble it up
    // (*numHeap)++ ;
    Heap->sizeHeap++ ;

    VSPINT nodeChild = Heap->sizeHeap, nodePar, valPar ;

    while(nodeChild > 1)
    {
        nodePar = nodeChild / 2 ;
        valPar  = heap[nodePar] ;
        if( valPar > valToInsert )
        {
            // Overwrite the child node with the parent node
            INSERT_IN_LIST(valPar, nodeChild, heap, ptrs) ;
            nodeChild = nodePar ;
        }
        else
        {
            break ;
        }
    }

    // Write i into to its proper location in the heap
    INSERT_IN_LIST(valToInsert, nodeChild, heap, ptrs) ;

    #ifdef VSP_DEBUG
        Min_Heap_Check_v4(Heap) ;
    #endif
    return ;
}



// ASSUMED: toDel is a value in the heap
// Remark: I think this function can be made more efficient by not actually
// writing lastVal until it finds its final home
void Min_Heap_Delete_v4
(
    Min_Heap_v4 *Heap,
    VSPINT      toDel // value to delete
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        Min_Heap_Check_v4(Heap) ;
    #endif

    VSPINT *heap  = Heap->heap ;
    VSPINT *ptrs  = Heap->heapPtrs ;

    // Overwrite toDel with the last value in the heap
    VSPINT heapNodeToDel    = ptrs[toDel] ;
    VSPINT lastVal          = heap[ (Heap->sizeHeap)-- ] ;
    heap[heapNodeToDel]     = lastVal ;
    ptrs[lastVal]           = heapNodeToDel ;
    ptrs[toDel]             = EMPTY ;

    // If we did NOT just delete the root node
    if(heapNodeToDel > 1 )
    {
        // Check if it needs to be up or down
        VSPINT parNode = heapNodeToDel / 2 ;
        VSPINT valPar = heap[parNode] ;
        if(lastVal < valPar )
        {
            // Bubble up
            Min_Heap_Bubble_Up_v4(Heap, heapNodeToDel, lastVal) ;
        }
        else if(heapNodeToDel <= Heap->sizeHeap)
        {
            // If we deleted the very last node, then heapNodeToDel would be
            // out of bounds for hea heap, which would cause a bug
            // Bubble down
            Min_Heap_Bubble_Down_v4(Heap, heapNodeToDel) ;
        }
    }
    else // We deleted the root node 
    {
        VSP_ASSERT(heapNodeToDel == 1) ;
        // Bubble down messes things up if we call it on an empty heap
        if(Heap->sizeHeap)
        {
            Min_Heap_Bubble_Down_v4(Heap, 1) ;
        }
    }

    #ifdef VSP_DEBUG
        Min_Heap_Check_v4(Heap) ;
    #endif
    return ;
}


/*----------------------------------------------------------------------------*/
/*   MAX HEAP FUNCTIONS                                                       */
/*----------------------------------------------------------------------------*/

Max_Heap_v4 *Max_Heap_Init_v4
(
    VSPINT gamma,
    VSPINT    *c,
    VSPINT     n  // number of vertices at current level
)
{
    Max_Heap_v4* heap = malloc(sizeof(*heap)) ;
    heap->c = c ;
    
    heap->heap        = malloc( ( gamma + 1 ) * sizeof(VSPINT)) ;
    heap->heapPtrs    = malloc( ( gamma + 1 ) * sizeof(VSPINT)) ;
    for(VSPINT i = 0 ; i < gamma + 1 ; i++)
    {
        heap->heapPtrs[i] = EMPTY ;
    }
    heap->sizeHeap    = 0 ;
    heap->gamma = gamma ;
    heap->n     = n ;

    return heap;
}

// Destructor for Max_Heap_v4
void Max_Heap_Destroy_v4(Max_Heap_v4 *heap)
{
    free(heap->heap);
    free(heap->heapPtrs);
    free(heap);
}

void Max_Heap_Check_v4
(
    Max_Heap_v4 *Heap
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT n = Heap->n ;
    VSPINT i, p, cval ;
    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;
    VSPINT  sizeHeap = Heap->sizeHeap ;
    VSPINT  gamma = Heap->gamma ;
    VSPINT *c = Heap->c ;
    VSP_ASSERT(n > 0) ;

    // Make sure size does not exceed maximum
    if (sizeHeap > gamma)
    {
        ERR("Heap size exceeds gamma") ;
    }

    /*------------------------------------------------------------------------*/
    // Check 
    // 1) Each value in the heap is in the correct range and it a c value
    // 2) The heap property is satisfies
    // 3) The pointer array points to it
    /*------------------------------------------------------------------------*/
    for(VSPINT par = 1 ; par <= sizeHeap ; par++)
    {
        VSPINT left   = 2 * par ;
        VSPINT right  = 2 * par + 1 ;
        VSPINT cpar   = heap[par] ;

        // Check that cpar is a value in the c array
        if( !Is_Array_Member(cpar, c, 0, n) )
        {
            ERR("Invalid entry into array") ;
        }
        
        // Check that the entry is in the valid range
        if( !IN_RANGE(cpar, 1, gamma + 1) )
        {
            ERR("Invalid entry in heap") ;
        }

        // Check the heap property
        if(left <= sizeHeap)
        {
            VSPINT cleft  = heap[left] ;
            if( !( cleft <= cpar ) )
            {
                ERR("Heap Property Violated") ;
            }
        }

        if(right <= sizeHeap)
        {
            VSPINT cright = heap[right] ;
            if( !( cright <= cpar ) )
            {
                ERR("Heap Property Violated") ;
            }
        }

        // Check the pointer
        if( ptrs[cpar] != par )
        {
            ERR("Heap pointers wrong") ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Check
    // 1) ptrs[cval] != 1 if and only if cval is not in the heap 
    /*------------------------------------------------------------------------*/
    for(VSPINT cval = 1 ; cval <= gamma ; cval++)
    {
        if(ptrs[cval] == EMPTY)
        {
            if( Is_Array_Member(cval, heap, 1, sizeHeap + 1) )
            {
                ERR("cval in heap but pointer says it's not") ;
            }
        }
        else
        {
            if( !Is_Array_Member(cval, heap, 1, sizeHeap + 1) )
            {
                ERR("Pointers say cval is in heap but it's actually not") ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Check the heap has no duplicated
    /*------------------------------------------------------------------------*/
    VSPINT *cvalsAppearing = calloc((gamma + 1), sizeof(VSPINT) ) ;
    for(p = 1 ; p <= sizeHeap ; p++)
    {
        cval = heap[p] ;
        cvalsAppearing[cval]++ ;
    }

    for(i = 0 ; i <= gamma ; i++)
    {
        if(cvalsAppearing[i] != 0 && cvalsAppearing[i] != 1)
        {
            ERR("Duplicate in heap") ;
        }
    }
    free(cvalsAppearing) ;
}



/* On input:
The subtree of "heap" starting at node with index "p" is a heap,
except that it's root may be out of order.

On output: 
subtree of "heap" starting at node with index "p" will be a subheap. */
void Max_Heap_Bubble_Down_v4
(
    Max_Heap_v4     *Heap,
    VSPINT              p  // start at node p in the heap
)
{
    VSPINT left, right, e, hleft, hright ;

    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;
    VSPINT nheap = Heap->sizeHeap ;

    e  = heap[p] ;
    while ( 1 )
    {
        left  = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            // Both left and right children of p are in the heap
            hleft  = heap[left] ;
            hright = heap[right] ;

            // We MUST swap with the larger child, there is no choice. If we
            // swap with the smaller child, then the new "root" will violate
            // the heap property because it will be smaller than one of its
            // children
            if ( hleft > hright )
            {
                // Traverse down the left heap
                if ( e < hleft )
                {
                    // Swap e with the left child
                    INSERT_IN_LIST(hleft, p, heap, ptrs) ;
                    p = left ;
                }
                else
                {
                    // Element e, with value ce, has found its final home
                    INSERT_IN_LIST(e, p, heap, ptrs) ;
                    return ;
                }
            }
            else
            {
                // Traverse down the right heap
                if ( e < hright )
                {
                    INSERT_IN_LIST(hright, p, heap, ptrs) ;
                    p = right ;
                }
                else
                {
                    // Element e, with value ce, has found its final home
                    INSERT_IN_LIST(e, p, heap, ptrs) ;
                    return ;
                }
            }
        }
        else
        {
            // Only the left child of p is in the heap
            if ( left <= nheap )
            {
                // Traverse down the left heap
                hleft = heap [left] ;
                if ( e < hleft )
                {
                    // Swap e with the left child
                    INSERT_IN_LIST(hleft, p, heap, ptrs) ;
                    p = left ;
                }
            }
            // Element e, with value ce, has found its final home
            INSERT_IN_LIST(e, p, heap, ptrs) ;
            return ;
        }
    }
}


// INPUT: The past lineage of p is all in order, i.e. 
// parent >= grandparent >= grand-grandparent .... 
// but p MAY be < parent
//
// Any children of p is ignored. They may or may not satisfy the heap property
//
// OUTPUT: the value residing at node p is in the correct place in its lineage
void Max_Heap_Bubble_Up_v4
(
    Max_Heap_v4         *Heap,
    VSPINT                  p, // node in the heap to bubble up
    VSPINT        valToBubble // value (equivalently, key) stored initially in node p
)
{
    VSPINT nodePar, valPar, keyPar ;
    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;

    while(p > 1)
    {
        nodePar = p / 2 ;
        valPar  = heap[nodePar] ;
        if( valPar < valToBubble )
        {
            // Overwrite the child node with the parent node
            INSERT_IN_LIST(valPar, p, heap, ptrs) ;
            p = nodePar ; // TODO: Consider p = p / 2
        }
        else
        {
            break ;
        }
    }

    // Write valToBubble into to its proper location in the heap
    INSERT_IN_LIST(valToBubble, p, heap, ptrs) ;
    return ;
}

// Takes the array "heap" and heapifies it
void Max_Heap_Build_v4
(
    Max_Heap_v4 *Heap
)
{
    VSPINT n = Heap->sizeHeap ;
    for (VSPINT p = n/2  ; p >= 1  ; p--)
    {
        Max_Heap_Bubble_Down_v4(Heap, p)  ;
    }
}

// Trick to make this more efficient: when we bubble up and swap two nodes,
// it is not neccesary to actually write the child node to the parent node.
// It is sufficient to write the parent node to the child node. At the very end,
// when i can no longer bubble up, it can be written to the current child
void Max_Heap_Insert_v4
(
    Max_Heap_v4         *Heap,
    VSPINT        valToInsert  // value/key to insert
)
{
    #ifdef VSP_DEBUG
        Max_Heap_Check_v4(Heap) ;
    #endif
    // Insert the value at the end and bubble it up
    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;

    if(ptrs[valToInsert] != EMPTY)
    {
        return ;
    }

    Heap->sizeHeap++ ;
    VSPINT nodeChild = Heap->sizeHeap, nodePar, valPar ;

    while(nodeChild > 1)
    {
        nodePar = nodeChild / 2 ;
        valPar  = heap[nodePar] ;
        if( valPar < valToInsert )
        {
            // Overwrite the child node with the parent node
            INSERT_IN_LIST(valPar, nodeChild, heap, ptrs) ;
            nodeChild = nodePar ;
        }
        else
        {
            break ;
        }
    }

    // Write i into to its proper location in the heap
    INSERT_IN_LIST(valToInsert, nodeChild, heap, ptrs) ;
    #ifdef VSP_DEBUG
        Max_Heap_Check_v4(Heap) ;
    #endif
    return ;
}


// ASSUMED: toDel is a value in the heap
// Remark: I think this function can be made more efficient by not actually
// writing lastVal until it finds its final home
void Max_Heap_Delete_v4
(
    Max_Heap_v4 *Heap,
    VSPINT      toDel  // value to delete
)
{
    SET_FUNCTION_COUNTER ;

    #ifdef VSP_DEBUG
        Max_Heap_Check_v4(Heap) ;
    #endif

    VSPINT *heap = Heap->heap ;
    VSPINT *ptrs = Heap->heapPtrs ;


    // Overwrite toDel with the last value in the heap
    VSPINT heapNodeToDel    = ptrs[toDel] ;
    VSPINT lastVal          = heap[ (Heap->sizeHeap)-- ] ;
    heap[heapNodeToDel]     = lastVal ;
    ptrs[lastVal]           = heapNodeToDel ;
    ptrs[toDel] = EMPTY ;

    // If we did NOT just delete the root node
    if(heapNodeToDel > 1 )
    {
        // Check if it needs to be up or down
        VSPINT parNode = heapNodeToDel / 2 ;
        VSPINT valPar = heap[parNode] ;
        if(lastVal > valPar )
        {
            // Bubble up
            Max_Heap_Bubble_Up_v4(Heap, heapNodeToDel, lastVal) ;
        }
        else if(heapNodeToDel <= Heap->sizeHeap) 
        {
            // If we deleted the very last node, then heapNodeToDel would be
            // out of bounds for hea heap, which would cause a bug
            // Bubble down
            Max_Heap_Bubble_Down_v4(Heap, heapNodeToDel) ;
        }
    }
    else // We deleted the root node 
    {
        VSP_ASSERT(heapNodeToDel == 1) ;
        // Bubble down messes things up if we call it on an empty heap
        if(Heap->sizeHeap)
        {
            Max_Heap_Bubble_Down_v4(Heap, 1) ;
        }
    }

    #ifdef VSP_DEBUG
        Max_Heap_Check_v4(Heap) ;
    #endif
    return ;
}

// Return top key in the heap. If the heap is empty, return EMPTY
VSPINT Max_Heap_See_Top_Key_v4
(
    Max_Heap_v4 *Heap
)
{
    if(Heap->sizeHeap)
    {
        return Heap->heap[1] ;
    }
    else
    {
        return EMPTY ;
    }
}
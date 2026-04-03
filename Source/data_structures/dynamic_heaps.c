#include "core.h"
#include "ds_api.h"

Dynamic_Heap *Dynamic_Heap_Init
(
    VSPINT                  *c,
    VSPINT               *ptrs,
    VSPINT          maxNumVals,
    VSPINT               minOrMax  // 1 for max, -1 for min
)
{
    SET_FUNCTION_COUNTER ;

    #ifdef VSP_DEBUG
        // Validate input
        TMP_ASSERT( c != NULL ) ;
        TMP_ASSERT( ptrs != NULL ) ;
        TMP_ASSERT( maxNumVals > 0 ) ;
        TMP_ASSERT( minOrMax == MIN_HEAP || minOrMax == MAX_HEAP ) ;
    #endif

    Dynamic_Heap *Heap = calloc( 1, sizeof(*Heap) ) ;
    Heap->heap = malloc( (maxNumVals + 1) * sizeof(VSPINT) ) ;
    Heap->ptrs = ptrs ;
    Heap->c = c ;
    Heap->numVals = 0 ;
    Heap->maxNumVals = maxNumVals ;
    Heap->minOrMax = minOrMax ;
    return Heap ;
}



void Dynamic_Heap_Destroy
(
    Dynamic_Heap *Heap
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        // Validate input
        TMP_ASSERT(Heap != NULL) ;
        TMP_ASSERT(Heap->heap != NULL) ;
    #endif


    free(Heap->heap) ;
    free(Heap) ;
    return ;
}


// Assumed: toInsert is NOT in the heap
void Dynamic_Heap_Insert
(
    Dynamic_Heap    *Heap,
    VSPINT       toInsert  // key to insert
)
{
    SET_FUNCTION_COUNTER ;

    #ifdef VSP_DEBUG
        // Check that toInsert is NOT in the heap
        if( Is_Array_Member(toInsert, Heap->heap, 1, Heap->numVals + 1) )
        {
            ERR("Trying to insert duplicate value in heap") ;
        }
    #endif

    // If the heap is full, allocate more memory
    if(Heap->numVals == Heap->maxNumVals)
    {
        // Don't forget the + 1
        VSPINT newMaxSize = 2 * Heap->maxNumVals ;
        Heap->heap = realloc( Heap->heap, (newMaxSize + 1) * sizeof(VSPINT) ) ;
        Heap->maxNumVals = 2 * Heap->maxNumVals ;
    }

    if(Heap->minOrMax == MAX_HEAP)
    {
        MaxHeap_Insert_v0
        (
            Heap->heap,
            Heap->ptrs,
            &(Heap->numVals),
            Heap->c,
            toInsert
        ) ;
    }
    else
    {
        MinHeap_Insert_v0
        (
            Heap->heap,
            Heap->ptrs,
            &(Heap->numVals),
            Heap->c,
            toInsert
        ) ;
    }
    return ;
}

// Returns EMPTY if heap is empty
void Dynamic_Heap_Get_Top_Val
(
    Dynamic_Heap   *Heap,
    VSPINT       *topVal
)
{
    if(Heap->numVals)
    {
        if(Heap->minOrMax == MAX_HEAP)
        {
            MaxHeap_Get_Top_Val_v0
            (
                Heap->heap,
                Heap->ptrs,
                &(Heap->numVals),
                Heap->c,
                topVal
            ) ;
        }
        else
        {
            MinHeap_Get_Top_Val_v0
            (
                Heap->heap,
                Heap->ptrs,
                &(Heap->numVals),
                Heap->c,
                topVal
            ) ;
        }

        // Free up some memory if the heap is empty
        if(Heap->numVals == 0)
        {
            // Don't forget the + 1
            Heap->heap = realloc(Heap->heap, (5 + 1) * sizeof(VSPINT) ) ;
            Heap->maxNumVals =  5 ;
        }
    }
    else
    {
        *topVal = EMPTY ;
    }
    return ;
}

// Returns EMPTY if heap is empty
void Dynamic_Heap_See_Top_Val
(
    Dynamic_Heap   *Heap,
    VSPINT       *topVal
)
{
    if(Heap->numVals)
    {
        *topVal = Heap->heap[1] ;
    }
    else
    {
        *topVal = EMPTY ;
    }
    return ;
}



// Need arbitrary deletion because of neighbors of vertex
// ASSUMED: toDel is in the heap
void Dynamic_Heap_Delete
(
    Dynamic_Heap *Heap,
    VSPINT       toDel
)
{
    if(Heap->minOrMax == MAX_HEAP)
    {
        MaxHeap_Delete_v0
        (
            Heap->heap,
            Heap->ptrs,
            &(Heap->numVals),
            Heap->c,
            toDel
        ) ;
    }
    else
    {
        MinHeap_Delete_v0
        (
            Heap->heap,
            Heap->ptrs,
            &(Heap->numVals),
            Heap->c,
            toDel
        ) ;
    }

    // Free up some memory if the heap is empty
    if( Heap->numVals == 0 )
    {
        // Don't forget the + 1
        Heap->heap = realloc(Heap->heap, (5 + 1) * sizeof(VSPINT)) ;
        Heap->maxNumVals = 5 ;
    }

    return ;
}

void Dynamic_Heap_Check
(
    Dynamic_Heap *Heap
)
{
    DBG_FN_BOILERPLATE ;
    TMP_ASSERT(Heap->maxNumVals > 0) ;
    

    if(Heap->minOrMax == MAX_HEAP)
    {
        MaxHeap_Check_v0(Heap->heap, Heap->ptrs, Heap->numVals, Heap->c) ;
    }
    else if(Heap->minOrMax == MIN_HEAP)
    {
        MinHeap_Check_v0(Heap->heap, Heap->ptrs, Heap->numVals, Heap->c) ;
    }
    else
    {
        ERR("Dynamic heap flag is neither min nor max") ;
    }

    if( !IN_RANGE(Heap->numVals, 0, Heap->maxNumVals + 1) )
    {
        ERR("Bad number of values in dynamic heap") ;
    }
    return ;
}
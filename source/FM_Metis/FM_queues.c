#include "FM_internal.h"



/*----------------------------------------------------------------------------*/
/*      Constructor/Destructor                                                */
/*----------------------------------------------------------------------------*/
QueueFM *QueueFM_Construct
(
    VSPINT maxNumNodes
)
{
    /* Construct the queue */
    QueueFM *q = malloc( sizeof(*q) ) ;
    q->maxNumNodes = maxNumNodes ;

    /* Initialize the heap */
    q->heap     = malloc( (maxNumNodes + 1) * sizeof(Node) ) ;
    q->numNodes = 0 ;

    /* Initialize the locator */
    q->locator = malloc( maxNumNodes * sizeof(VSPINT) ) ;
    // TODO: Consider if placing q->locator in a local var speeds things up
    for(VSPINT i = 0 ; i < maxNumNodes ; i++) q->locator[i] = EMPTY ;

    return q ;
}

void QueueFM_Destruct
(
    QueueFM *q
)
{
    free(q->heap) ;
    free(q->locator) ;
    free(q) ;
    return ;
}

/*----------------------------------------------------------------------------*/
/*      Top Key/Value                                                         */
/*----------------------------------------------------------------------------*/
VSPINT QueueFM_See_Top_Value(QueueFM *q)
{
    if(q->numNodes == 0) return EMPTY ;
    else return q->heap[1].value ;
}

/* Assumes queue is NOT empty */
VSPINT QueueFM_See_Top_Key(QueueFM *q)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(q->numNodes < 1)
        {
            VSP_error("Acessing top key of an empty FM queue", DBG_INFO) ;
        }
    #endif

    return q->heap[1].key ;
}

VSPINT QueueFM_Get_Top_Value(QueueFM *q)
{
    VSPINT top = QueueFM_See_Top_Value(q) ;
    if(top != EMPTY)
    {
        QueueFM_Delete(q, top) ;
    }

    return top ;
}

/*----------------------------------------------------------------------------*/
/*      Insert/Delete                                                         */
/*----------------------------------------------------------------------------*/
/* Works as follows: create an empty node at the end of the heap. Examine the
   parent node, and if its key is smaller than the new key, push the parent
   down to fill the hole. Repeat until either the key of the parent of the hole
   is g.e. the new key (at which point you place the new value and key in the
   hole), or until the hole is at the root, at which point you also place the
   new key and value into the hole 
   
   ASSUMED: newValue is NOT a member of the heap already
   */
void QueueFM_Insert(QueueFM *q, VSPINT newValue, VSPINT newKey)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(q->locator[newValue] != EMPTY)
        {
            VSP_error("Trying to insert duplicate value in FM heap", DBG_INFO) ;
        }

        QueueFM_Check_Heap_Light(q) ;
    #endif
    /* Unpack queue into local vars */
    Node *heap       = q->heap ;
    VSPINT *locator  = q->locator ;

    /* Initalize child and parent nodes */
    VSPINT childNode = (++q->numNodes) ;

    /* Bubble the new node up */
    while( childNode > 1 )
    {
        VSPINT parNode = childNode / 2 ;
        if(heap[parNode].key < newKey )
        {
            heap[childNode] = heap[parNode] ;
            locator[ heap[childNode].value ] = childNode ;

            childNode = parNode ;
        }
        else break ;
    }
    heap[childNode].value = newValue ;
    heap[childNode].key   = newKey ;
    locator[newValue]     = childNode ;

    #ifdef VSP_DEBUG
        QueueFM_Check_Heap_Light(q) ;
    #endif

    return ;
}

/* ASSUME: The value is present in the heap */
void QueueFM_Delete(QueueFM *q, VSPINT valueToDelete)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        QueueFM_Check_Heap_Light(q) ;
    #endif

    /* Unpack queue into local variables */
    Node   *heap     = q->heap;
    VSPINT *locator  = q->locator;

    /* Update the locator */
    VSPINT oldNode         = locator[valueToDelete] ;
    locator[valueToDelete] = EMPTY ;
    #ifdef VSP_DEBUG
        /* Check that the value is in the heap */
        if(oldNode == EMPTY)
        {
            VSP_error("Trying to delete node which is not in the queue", DBG_INFO) ;
        }
    #endif

    /* Update numNodes and quit early if the value to delete was in the very
       last node. Note that this early exit also catches the case when the
       heap is now empty */
    if( heap[ q->numNodes-- ].value == valueToDelete )
    {
        #ifdef VSP_DEBUG
            QueueFM_Check_Heap_Light(q) ;
        #endif
        return ;
    }

    VSPINT  i, j, numNodes ;
    VSPINT newValue = heap[ q->numNodes + 1 ].value ; 
    VSPINT newKey   = heap[ q->numNodes + 1 ].key ; 
    VSPINT oldKey   = heap[oldNode].key ;

    i        = oldNode ;
    numNodes = q->numNodes ;

    if(newKey < oldKey) // Bubble down if necessary
    {
        while( (j = 2 * i) <= numNodes)
        {
            if(heap[j].key > newKey)
            {
                /* Check if right child is EVEN GREATER than left child */
                if(j + 1 <= numNodes && heap[j + 1].key > heap[j].key) j++ ;
                heap[i] = heap[j] ;
                locator[heap[i].value] = i;
                i = j;
            }
            else if(j + 1 <= numNodes && heap[j + 1].key > newKey)
            {
                j++ ;
                heap[i] = heap[j] ;
                locator[heap[i].value] = i;
                i = j;
            }
            else break ;
        }
    }
    else // Bubble up
    {
        while( i > 1 )
        {
            j = i / 2 ;

            /* If the parent is smaller than newKey, swap */
            if( heap[j].key < newKey )
            {
                heap[i] = heap[j] ;
                locator[heap[i].value] = i ;
                i = j;
            }
            else break ;
        }
    }

    heap[i].value     = newValue ;
    heap[i].key       = newKey ;
    locator[newValue] = i ;

    #ifdef VSP_DEBUG
        QueueFM_Check_Heap_Light(q) ;
    #endif

    return;
}

/*----------------------------------------------------------------------------*/
/*      Update                                                                */
/*----------------------------------------------------------------------------*/
/* ASSUMED: value is in the heap */
void QueueFM_Update(QueueFM *q, VSPINT value, VSPINT newKey)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        if(value < 0 || value > q->maxNumNodes - 1)
        {
            VSP_error("Invalid node", DBG_INFO) ;
        }

        if(q->locator[value] == EMPTY)
        {
            VSP_error("Value not in heap", DBG_INFO) ;
        }
        QueueFM_Check_Heap_Light(q) ;
    #endif

    /* Unpack the struct into local variables */
    Node   *heap        = q->heap ;
    VSPINT *locator     = q->locator ;
    VSPINT  numNodes    = q->numNodes ;

    VSPINT i, j ;
    i = locator[value] ;
    VSPINT oldKey = heap[i].key ;

    if ( newKey > oldKey )
    { 
        /* Bubble up */
        while (i > 1)
        {
            j = i / 2 ;

            /* If new key is bigger than parent key, swap */
            if( newKey > heap[j].key )
            {
                heap[i] = heap[j];
                locator[heap[i].value] = i;
                i = j;
            }
            else break ;
        }
    }
    else
    { 
        /* Bubble down */
        while ( (j = 2 * i) <= numNodes )
        {
            /* If the new is smaller than child, swap */
            if ( newKey < heap[j].key)
            {
                /* Swap with the larger child */
                if (j + 1 <= numNodes && heap[j + 1].key > heap[j].key ) j++ ;
                heap[i] = heap[j];
                locator[heap[i].value] = i;
                i = j;
            }
            else if (j + 1 <= numNodes && newKey < heap[j + 1].key)
            {
                j++;
                heap[i] = heap[j];
                locator[heap[i].value] = i;
                i = j;
            }
            else
                break;
        }
    }

    heap[i].key = newKey;
    heap[i].value = value;
    locator[value] = i;

    #ifdef VSP_DEBUG
        QueueFM_Check_Heap_Light(q) ;
    #endif
}




/*----------------------------------------------------------------------------*/
/*      Debugging                                                             */
/*----------------------------------------------------------------------------*/
/* Check as much as possible while only having access to the data in the queue*/
void QueueFM_Check_Heap_Light(QueueFM *q)
{
    static long long unsigned count = 0 ;
    count++ ;

    /* Unpack the struct into local variables */
    Node   *heap        = q->heap ;
    VSPINT *locator     = q->locator ;
    VSPINT  numNodes    = q->numNodes ;
    VSPINT  maxNumNodes = q->maxNumNodes ;

    /* Check that numNodes is non-negative */
    if(q->numNodes < 0)
    {
        char msg[ERR_MSG_LEN] ;
        sprintf(msg, "Num nodes in FM heap: %i, which is negative",
                     q->numNodes) ;
        VSP_error(msg, DBG_INFO) ;
    }

    /* Check that maxNumNodes is positive */
    if(q->maxNumNodes <= 0)
    {
        char msg[ERR_MSG_LEN] ;
        sprintf(msg, "Max num nodes in FM heap: %i, which is non-positive",
                     q->maxNumNodes) ;
        VSP_error(msg, DBG_INFO) ;
    }

    /* Check that numNodes does not exceed maxNumNodes */
    if(q->numNodes > q->maxNumNodes)
    {
        char msg[ERR_MSG_LEN] ;
        sprintf(msg, "Num nodes in FM heap: %i, which is greater than the"
                     "maximum, %i", q->numNodes, q->maxNumNodes              ) ;
        VSP_error(msg, DBG_INFO) ;
    }

    for(VSPINT i = 1 ; i <= q->numNodes ; i++)
    {
        // Check that the value is in the correct range
        if( !(0 <= heap[i].value && heap[i].value < q->maxNumNodes) )
        {
            VSP_error("Out of bounds value in the heap", DBG_INFO) ;
            break ;
        }

        /* Check that left child is smaller */
        VSPINT child = 2 * i ;
        if(child <= q->numNodes)
        {
            if(heap[i].key < heap[child].key)
            {
                VSP_error("Heap property violated in an FM queue", DBG_INFO) ;
                break ;
            }
        }

        /* Check that right child is smaller */
        child++ ;
        if(child <= q->numNodes)
        {
            if(heap[i].key < heap[child].key)
            {
                VSP_error("Heap property violated in an FM queue", DBG_INFO) ;
                break ;
            }
        }

        /* Check that the locator is correct */
        if(locator[heap[i].value] != i)
        {
            VSP_error("Locator is incorrect in FM queue", DBG_INFO) ;
            break ;
        }
    }

    /* Check that the locator does not contain too many values */
    for(VSPINT i = 0 ; i < maxNumNodes ; i++)
    {
        /* Determine if i is in the heap */
        VSPINT k ;
        bool is_i_in_heap = false ;
        for(k = 1 ; k <= numNodes ; k++)
        {
            if(heap[k].value == i)
            {
                is_i_in_heap = true ;
                break ;
            }
        }

        if(is_i_in_heap)
        {
            if(locator[i] != k)
            {
                VSP_error("Locator in FM queue is incorrect", DBG_INFO) ;
            }
        }
        else
        {
            if(locator[i] != EMPTY)
            {
                VSP_error("Locator in FM queue says a node is in"
                         "the heap, when it's not", DBG_INFO) ;
            }
        }
    }

    /* Check that there are no duplicate values */
    VSPINT *valsSeen = calloc((q->maxNumNodes), sizeof(VSPINT) ) ;
    for(VSPINT i = 1 ; i <= q->numNodes ; i++)
    {
        valsSeen[i]++ ;
        if(valsSeen[i] != 0 && valsSeen[1] != 1)
        {
            VSP_error("Duplicate value in heap", DBG_INFO) ;
        }
    }
    free(valsSeen) ;
}
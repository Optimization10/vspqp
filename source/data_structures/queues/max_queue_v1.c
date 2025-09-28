#include "shared_api.h"
#include "ds_api.h"


static void VSP_max_heapify_v1
(
    Max_Queue_v1 *queuestr,
    VSPINT               p  /* start at node p in the heap */
) ;

VSPINT Max_Queue_See_Top_Entry_v1
(
    Max_Queue_v1 *queuestr
)
{
    if(queuestr->n_queue == 0)
    {
        return EMPTY ;
    }
    return queuestr->queue[1] ;
}

VSPINT VSP_max_queue_see_top_priority_v1
(
    Max_Queue_v1 *queuestr
)
{
    if(queuestr->n_queue == 0)
    {
        return EMPTY ;
    }
    return queuestr->priorities[ queuestr->queue[1] ] ;
}


/* Construct an uninitialized min priority queue v1 ds. Everything including the
   lookup is initialized */
Max_Queue_v1 *VSP_construct_max_queue_v1
(
    VSPINT n_priorities,
    VSPINT *priorities
)
{
    Max_Queue_v1 *ret = malloc(sizeof(*ret)) ;

    ret->queue        = malloc((n_priorities + 1) * sizeof(*(ret->queue))) ;
    ret->priorities   = priorities ;
    ret->n_queue      = 0 ;
    ret->n_priorities = n_priorities ;
    ret->lookup       = malloc((n_priorities) * sizeof(*(ret->lookup))) ;
    for(VSPINT i = 0 ; i < n_priorities ; i++)
    {
        ret->lookup[i] = EMPTY ;
    }
    return ret ;
}

VSPINT Max_Queue_See_Top_Key_v1
(
    Max_Queue_v1 *queue
)
{
    if(queue->n_queue == 0)
    {
        return EMPTY ;
    }
    else
    {
        return queue->priorities[ queue->queue[1] ] ;
    }
}


/* Destroy everything EXCEPT the priorities array */
void VSP_destruct_max_queue_v1
(
    Max_Queue_v1 *to_destroy
)
{
    free(to_destroy->queue) ;
    free(to_destroy->lookup) ;
    free(to_destroy) ;
    return ;
}


/* This function heapifies the queue array. It assumes the lookup array is 
   correct */
void VSP_max_queue_build_v1
(
    Max_Queue_v1 *queuestr
)
{
    VSPINT n_queue = queuestr->n_queue ;
    VSPINT p  ;

    for (p = n_queue/2  ; p >= 1  ; p--)
    {
        VSP_max_heapify_v1(queuestr, p) ;
    }

    #ifdef VSP_DEBUG
        // VSP_check_min_max_heap_v1(queuestr, MAX_HEAP) ;
    #endif

    return ;
}

// Remove the value "rm" from the queue
void Max_Queue_Delete_v1
(
    Max_Queue_v1 *queuestr,
    VSPINT rm
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        // VSP_check_min_max_heap_v1(queuestr, MAX_HEAP) ;
    #endif

    VSPINT *lookup = queuestr->lookup ;
    VSPINT address  = lookup[rm] ;
    if(address == EMPTY)
    {
        return ;
    }
    else if(address == queuestr->n_queue) /* rm is in the last node */
    {
        queuestr->n_queue-- ; ;
        lookup[rm] = EMPTY ;
        return ;
    }

    VSPINT *queue      = queuestr->queue ;
    VSPINT *priorities = queuestr->priorities ;
    VSPINT replacement = queue[ (queuestr->n_queue)-- ] ;
    lookup[rm] = EMPTY ;

    /* Determine if we need to bubble the replacement up or down */

    if(address == 1 || priorities[ queue[address / 2] ] >= priorities[replacement])
    {
        /* Bubble down */
        queue[address] = replacement ;
        lookup[replacement] = address ;
        VSP_max_heapify_v1(queuestr, address) ;
    }
    else
    {
        /* Bubble up */
        while(address > 1)
        {
            VSPINT par = address / 2 ;
            if(priorities[ queue[par] ] < priorities[replacement])
            {
                queue[address] = queue[par] ;
                lookup[queue[par]] = address ;
                address = par ;
            }
            else
            {
                break ;
            }
        }
        queue[address] = replacement ;
        lookup[replacement] = address ;
    }

    #ifdef VSP_DEBUG
        // VSP_check_min_max_heap_v1(queuestr, MAX_HEAP) ;
    #endif

    return ;
}

/* This function does not check that the element to add is not already in
   the queue */
void VSP_max_queue_add_v1
(
    Max_Queue_v1 *queuestr,
    VSPINT toadd
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT par, priorpar ;
    VSPINT *priorities = queuestr->priorities ;
    VSPINT *heap       = queuestr->queue ;
    VSPINT *lookup     = queuestr->lookup ;
    VSPINT child       = ++(queuestr->n_queue)  ;
    VSPINT priorchild  = priorities[toadd] ;

    #ifdef VSP_DEBUG
        if(lookup[toadd] != EMPTY)
        {
            VSP_error("Element already in queue", DBG_INFO) ;
        }
    #endif

    /* toadd is perpetually considered the child, until its final place
       is found */
    while(child > 1)
    {
        par = child / 2 ;
        priorpar = priorities[ heap[par] ] ;
        if(priorchild > priorpar)
        {
            /* Write the parent to the child */
            heap[child] = heap[par] ;
            lookup[ heap[par] ] = child ;
            child = par ;
        }
        else
        {
            break ;
        }
    }

    /* Write the child */
    heap[child] = toadd ;
    lookup[toadd] = child ;
    return ;
}



/* On input:
The subtree of the heap starting at node with index "p" is a heap,
except that it's root may be out of order.

On output: 
subtree of the heap starting at node with index "p" will be a subheap. */
static void VSP_max_heapify_v1
(
    Max_Queue_v1 *queuestr,
    VSPINT            p          /* start at node p in the heap */
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT *heap   = queuestr->queue ;
    VSPINT *X      = queuestr->priorities ;
    VSPINT *lookup = queuestr->lookup ;
    VSPINT  nheap  = queuestr->n_queue ;

    VSPINT Xe, Xleft, Xright ;
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
            if ((Xleft > Xright)) /* for min heap: (Xleft < Xright) */
            {
                /* traverse down the left heap */
                if ((Xe < Xleft)) /* for min heap: (Xe > Xleft) */
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
            if ( left <= nheap )
            {
                /* traverse down the left heap */
                hleft = heap [left] ;
                Xleft = X [hleft] ;
                if ((Xe < Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    lookup [hleft] = p ;
                    p = left ;
                }
            }
            /* element e, with value Xe, has found its final home */
            heap [p]   = e ;
            lookup [e] = p ;
            return ;
        }
    }
}

void Max_Queue_Reset_v1
(
    Max_Queue_v1 *qstr
)
{
    VSPINT n_queue = qstr->n_queue ;
    for(VSPINT i = 1 ; i <= n_queue ; i++)
    {
        VSPINT j = qstr->queue[i] ;
        qstr->lookup[j] = EMPTY ;
    }
    qstr->n_queue = 0 ;
}



// ASSUMED: to_insert is NOT in the queue
void Max_Queue_Insert_v1
(
    Max_Queue_v1      *qstr,
    VSPINT        to_insert
)
{
    SET_FUNCTION_COUNTER ;


    VSPINT node_par ;
    VSPINT *heap       = qstr->queue ;
    VSPINT *lookup     = qstr->lookup ;
    VSPINT node_child       = ++(qstr->n_queue)  ;
    VSPINT  *priorities = qstr->priorities ;
    VSPINT priorchild  = priorities[to_insert] ;

    #ifdef VSP_DEBUG
        if(lookup[to_insert] != EMPTY)
        {
            VSP_error("Element already in queue", DBG_INFO) ;
        }
    #endif

    /* to_insert is perpetually considered the child, until its final place
       is found */
    while(node_child > 1)
    {
        node_par = node_child / 2 ;
        if( priorities[ heap[node_par] ] < priorchild )
        {
            /* Write the parent to the child */
            heap[node_child] = heap[node_par] ;
            lookup[ heap[node_par] ] = node_child ;
            node_child = node_par ;
        }
        else
        {
            break ;
        }
    }

    /* Write the child */
    heap[node_child] = to_insert ;
    lookup[to_insert] = node_child ;

    return ;
}


void Max_Queue_Print_Node_v1
(
    Max_Queue_v1   *q,
    VSPINT       node,
    char         *msg  // Ignored if NULL
)
{
    int i ;
    FILE *file = fopen("output/log_debug.txt", "a") ;
    fprintf(file, "\n    Printing Heap Node:") ;
    if(msg != NULL)
    {
        fprintf(file, " %s", msg) ;
    }
    fprintf(file, "\n") ;

    if(node > q->n_queue)
    {
        fprintf(file, "Node out of bounds") ;
        fclose(file) ;
        return ;
    }


    // Get vtx and ci of the node
    VSPINT vtx = q->queue[node] ;
    VSPINT ci = q->priorities[vtx] ;

    // Figure out the number of minuses
    // The .01 is in case the argument is a power of 10.
    int num_digits_node = (int) ceil( log(node + 0.1) / log(10) ) ;
    int num_digits_vtx  = (int) ceil( log(vtx + 0.1) / log(10) ) ;
    int num_digits_ci   = (int) ceil( log(ci + 0.1) / log(10) ) ;
    int num_minuses = 19 + 4 + num_digits_node + num_digits_vtx + num_digits_ci ;

    // Print top row of box
    fprintf(file, "+") ;
    for(i = 0 ; i < num_minuses ; i++) fprintf(file, "-") ;
    fprintf(file, "+") ;

    // Print middle row of box
    fprintf(file, "\n|  node: %i  vtx: %i  ci: %i  |\n", node, vtx, ci) ;

    //Print bottom row of box
    fprintf(file, "+") ;
    for(i = 0 ; i < num_minuses ; i++) fprintf(file, "-") ;
    fprintf(file, "+") ;

    fclose(file) ;
}
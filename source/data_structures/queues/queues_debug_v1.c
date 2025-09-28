#include "shared_api.h"
#include "ds_api.h"

void VSP_check_min_heap_v1
(
    Min_Queue_v1* queuestr
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT *queue = queuestr->queue ;
    VSPINT *priorities = queuestr->priorities ;
    VSPINT *lookup = queuestr->lookup ;
    VSPINT n_queue = queuestr->n_queue ; 
    VSPINT i, p, left, right, vleft, vright, vpar ;

    /* Check heap property */
    for(p = 1 ; p <= n_queue ; p++)
    {
        left  = 2 * p ;
        right = 2 * p + 1 ;

        if(left > n_queue && right > n_queue)
        {
            break ;
        }

        vpar  = priorities[queue[p]] ;
        if(left <= n_queue)
        {
            vleft = priorities[queue[left]] ;
            if(vpar > vleft)
            {
                VSP_error("Heap property violated", DBG_INFO) ;
                return ;
            }
        }

        if(right <= n_queue)
        {
            vright = priorities[queue[right]] ;
            if(priorities[queue[p]] > vright)
            {
                VSP_error("Heap property violated", DBG_INFO) ;
                return ;
            }
        }
    }

    /* Check lookup array */
    VSPINT n_prior = queuestr->n_priorities ;
    VSPINT *dbg_lookup = malloc(n_prior * sizeof(*dbg_lookup)) ;

    /* Initialize dbg_lookup to all EMPTY's */
    for(i = 0 ; i < n_prior ; i++)
    {
        dbg_lookup[i] = EMPTY ;
    }

    /* Build dbg_lookup */
    for(p = 1 ; p <= n_queue ; p++)
    {
        dbg_lookup[queue[p]] = p ;
    }

    /* Check dbg_lookup against lookup */
    for(i = 0 ; i < n_prior ; i++)
    {
        if(lookup[i] != dbg_lookup[i])
        {
            VSP_error("Lookup array incorrect", DBG_INFO) ;
            break ;
        }
    }

    free(dbg_lookup) ;
    return ;
}


void VSP_check_max_heap_v1
(
    Max_Queue_v1* queuestr
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT *queue = queuestr->queue ;
    VSPINT *priorities = queuestr->priorities ;
    VSPINT *lookup = queuestr->lookup ;
    VSPINT n_queue = queuestr->n_queue ; 
    VSPINT i, p, left, right, vleft, vright, vpar ;

    /* Check heap property */
    for(p = 1 ; p <= n_queue ; p++)
    {
        left  = 2 * p ;
        right = 2 * p + 1 ;


        if(left > n_queue && right > n_queue)
        {
            break ;
        }

        vpar  = priorities[queue[p]] ;
        if(left <= n_queue)
        {
            vleft = priorities[queue[left]] ;
            if(vpar < vleft)
            {
                VSP_error("Heap property violated", DBG_INFO) ;
                return ;
            }
        }

        if(right <= n_queue)
        {
            vright = priorities[queue[right]] ;
            if(priorities[queue[p]] < vright)
            {
                VSP_error("Heap property violated", DBG_INFO) ;
                return ;
            }
        }
    }

    /* Check lookup array */
    VSPINT n_prior = queuestr->n_priorities ;
    VSPINT *dbg_lookup = malloc(n_prior * sizeof(*dbg_lookup)) ;

    /* Initialize dbg_lookup to all EMPTY's */
    for(i = 0 ; i < n_prior ; i++)
    {
        dbg_lookup[i] = EMPTY ;
    }

    /* Build dbg_lookup */
    for(p = 1 ; p <= n_queue ; p++)
    {
        dbg_lookup[queue[p]] = p ;
    }

    /* Check dbg_lookup against lookup */
    for(i = 0 ; i < n_prior ; i++)
    {
        if(lookup[i] != dbg_lookup[i])
        {
            VSP_error("Lookup array incorrect", DBG_INFO) ;
            break ;
        }
    }

    free(dbg_lookup) ;
    return ;
}
#include "core.h"
#include "shared_api.h"
#include "ds_api.h" 

/*----------------------------------------------------------------------------*/
/* Definitions */
/*----------------------------------------------------------------------------*/

/* Construct a dynamic stack. Initial size is 32 elements */
VSP_dyn_stk *VSP_construct_dyn_stk()
{
    VSP_dyn_stk *ds  = malloc(sizeof(*ds)) ;
    ds->max_size = 32 ;
    ds->stack     = malloc( (ds->max_size) * sizeof(VSPINT) ) ;
    ds->n_stack   = 0 ;

    return ds ;
}

void VSP_destruct_dyn_stk(VSP_dyn_stk* ds)
{
    free(ds->stack) ;
    free(ds) ;
}

void VSP_dyn_stk_push(VSP_dyn_stk* ds, VSPINT to_push)
{
    if (ds->n_stack == ds->max_size)
    {
        ds->max_size *= 2 ;
        ds->stack = realloc(ds->stack, ds->max_size * sizeof(VSPINT)) ;
    }

    ds->stack[ds->n_stack] = to_push ;
    ds->n_stack++ ;
}

bool VSP_dyn_stk_is_member(VSPINT i, VSP_dyn_stk* ds)
{
    for (VSPINT j = 0 ; j < ds->n_stack ; j++)
    {
        if (ds->stack[j] == i)
        {
            return true ;
        }
    }

    return false ;
}

/* Create a deep copy */
VSP_dyn_stk *VSP_dyn_stk_copy(VSP_dyn_stk *stack_in)
{
    VSP_dyn_stk *stack_out = VSP_construct_dyn_stk() ;

    for (VSPINT i = 0 ; i < stack_in->n_stack ; i++)
    {
        VSP_dyn_stk_push(stack_out, stack_in->stack[i]) ;
    }

    return stack_out ;
}

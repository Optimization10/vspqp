#include "core.h"
#include "options.h"
#include "VSPQP.h"
#include "stdio.h"




VSP_levels* VSP_construct_levels
(
    VSP_levels_selection_type selection_type
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSP_levels *levels = malloc(sizeof(*levels)) ;
    levels->top_levels                = false ;
    levels->bot_levels                = false ;
    levels->all_levels                = false ;
    levels->no_levels                 = false ;
    levels->specify_selected_levels   = false ;
    levels->specify_unselected_levels = false ;

    levels->levels_selected       = NULL ;
    levels->levels_not_selected   = NULL ;

    if(selection_type == TOP_LEVELS)
    {
        levels->top_levels     = true ;
        levels->num_top_levels = 0 ;
    }
    else if(selection_type == BOT_LEVELS)
    {
        levels->bot_levels     = true ;
        levels->num_bot_levels = 0 ;
    }
    else if(selection_type == ALL_LEVELS)
    {
        levels->all_levels = true ;
    }
    else if(selection_type == NO_LEVELS)
    {
        levels->no_levels = true ;
    }
    else if(selection_type == SPECIFY_SELECTED_LEVELS)
    {
        levels->specify_selected_levels = true ;
        levels->levels_selected = VSP_construct_dyn_stk() ;
    }
    else if(selection_type == SPECIFY_UNSELECTED_LEVELS)
    {
        levels->specify_unselected_levels = true ;
        levels->levels_not_selected = VSP_construct_dyn_stk() ;
    }
    else
    {
        #ifdef VSP_DEBUG
            Error("Invalid selection type", DBG_INFO) ;
        #endif
    }

    return levels ;
}


void VSP_destruct_levels
(
    VSP_levels* levels
)
{
    if(levels->specify_selected_levels)
    {
        VSP_destruct_dyn_stk(levels->levels_selected) ;
    }
    else if(levels->specify_unselected_levels)
    {
        VSP_destruct_dyn_stk(levels->levels_not_selected) ;
    }

    free(levels) ;
}


bool VSP_levels_is_member
(
    Problem *problem,
    VSP_levels  *levels
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;

        /* Check that only one  value is true */
        VSPINT sum =    levels->top_levels +
                        levels->bot_levels +
                        levels->all_levels +
                        levels->no_levels  +
                        levels->specify_selected_levels +
                        levels->specify_unselected_levels ;
        if(sum != 1)
        {
            Error("Invalid selection type", DBG_INFO) ;
        }
    #endif
    VSPINT problem_level = problem->level ;

    if(levels->all_levels)
    {
        return true ;
    }
    else if(levels->no_levels)
    {
        return false ;
    }
    else if(levels->top_levels)
    {
        return (problem_level < levels->num_top_levels) ;
    }
    else if(levels->bot_levels)
    {
        /* Find the bottom level */
        VSPINT bot_level = problem->level ;
        Problem *p = problem ;
        while(p->child_problem != NULL)
        {
            p = p->child_problem ;
            bot_level = p->level ;
        }
        return (problem_level > bot_level - levels->num_bot_levels) ;
    }
    else if(levels->specify_selected_levels)
    {
        bool ret = VSP_dyn_stk_is_member
                    (
                        problem_level,
                        levels->levels_selected
                    ) ;
        return ret ;
    }
    else if(levels->specify_unselected_levels)
    {
        bool ret = !VSP_dyn_stk_is_member
                        (
                            problem_level,
                            levels->levels_not_selected
                        ) ;
        return ret ;
    }
    else
    {
        #ifdef VSP_DEBUG
            Error("Invalid selection type", DBG_INFO) ;
        #endif
    }
    return false ;
}

void VSP_levels_select
(
    VSPINT       level_num,
    VSP_levels  *levels
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if( levels->specify_selected_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef VSP_DEBUG
        if(VSP_dyn_stk_is_member(level_num, levels->levels_selected))
        {
            Error("Level already selected", DBG_INFO) ;
        }
    #endif

    VSP_dyn_stk_push(levels->levels_selected, level_num) ;
}


void VSP_levels_unselect
(
    VSPINT      level_num,
    VSP_levels    *levels
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if( levels->specify_unselected_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef VSP_DEBUG
        if( VSP_dyn_stk_is_member
            (
                level_num,
                levels->levels_not_selected
            ) )
        {
            Error("Level already unselected", DBG_INFO) ;
        }
    #endif

    VSP_dyn_stk_push(levels->levels_not_selected, level_num) ;

    return ;
}


void VSP_levels_select_top
(
    VSPINT n_top_levels,
    VSP_levels *levels
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(levels->top_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    levels->num_top_levels = n_top_levels;

    return ;
}

void VSP_levels_select_bot
(
    VSPINT n_bot_levels,
    VSP_levels *levels
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        if(levels->bot_levels == false)
        {
            Error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    levels->num_bot_levels = n_bot_levels;

    return ;
}




/*----------------------------------------------------------------------------*/
/* Dynamic Stack Definitions */
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





/* Just construct and set everything to empty */
VSP_restart_decision *VSP_construct_restart_decision()
{
    VSP_restart_decision *ret = malloc(sizeof(*ret)) ;

    /* Levels */
    ret->use_levels = false ;
    ret->levels     = NULL ;

    /* Rule */
    ret->use_rule             = false ;
    ret->restarts_done        = 0 ;
    ret->restarts_allowed     = 0 ;
    ret->rule_param           = 0 ;
    ret->finest_level_allowed = 0 ;

    return ret ;
}

void VSP_destruct_restart_decision(VSP_restart_decision *dec)
{
    if (dec->levels != NULL)
    {
        VSP_destruct_levels(dec->levels) ;
    }

    free(dec) ;

    return ;
}

/* If using the rule, and if it returns true, then it also increments
   dec->restarts_done. */
bool VSP_restart_decide
(
    VSP_restart_decision       *dec,
    VSPFLOAT                  objIn, /* If not using rule, this is ignored */
    VSPFLOAT                 objOut, /* If not using rule, this is ignored */
    Problem            *problem  /* If using levels, this is ignored   */
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        if(dec->use_levels + dec->use_rule != 1)
        {
            Error("Error in restart decision", DBG_INFO) ;
        }

        if(dec->use_rule)
        {
            if(dec->restarts_done > dec->restarts_allowed)
            {
                Error("Error in restart decision", DBG_INFO) ;
            }

            if(objOut < objIn)
            {
                Error("Probable error. Final objective worse.", DBG_INFO) ;
            }
        }
        else
        {
            if(dec->levels == NULL)
            {
                Error("Error in restart decision", DBG_INFO) ;
            }
        }
    #endif



    if (dec->use_levels)
    {
        return VSP_levels_is_member(problem, dec->levels) ;
    }
    else if (dec->use_rule)
    {
        if(dec->restarts_done >= dec->restarts_allowed
            || problem->level < dec->finest_level_allowed)
        {
            return false ;
        }

        /* Check if we should restart */
        VSPFLOAT rel_change = (objOut - objIn) / objIn ;
        if (rel_change <= dec->rule_param)
        {
            dec->restarts_done++ ;
            return true ;
        }
        else
        {
            return false ;
        }
    }

    return false ;
}

void VSP_restart_decide_reset(VSP_restart_decision *dec)
{
    dec->restarts_done = 0 ;
    return ;
}



/***************** Options ************************************************/
Options* VSP_construct_options()
{
    Options* options = malloc(sizeof(Options));
    options->restart_decision     = NULL ;
    return options;
}


void VSP_destruct_options(Options* options)
{
    if(options->restart_decision != NULL)
    {
        VSP_destruct_restart_decision(options->restart_decision) ;
    }


    /* Free the struct itself */
    free(options);
}
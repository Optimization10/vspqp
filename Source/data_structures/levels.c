#include "core.h"
#include "shared_api.h"
#include "ds_api.h" // IWYU pragma: keep


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
            VSP_error("Invalid selection type", DBG_INFO) ;
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
    VSP_problem *problem,
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
            VSP_error("Invalid selection type", DBG_INFO) ;
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
        VSP_problem *p = problem ;
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
            VSP_error("Invalid selection type", DBG_INFO) ;
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
            VSP_error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef VSP_DEBUG
        if(VSP_dyn_stk_is_member(level_num, levels->levels_selected))
        {
            VSP_error("Level already selected", DBG_INFO) ;
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
            VSP_error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    #ifdef VSP_DEBUG
        if( VSP_dyn_stk_is_member
            (
                level_num,
                levels->levels_not_selected
            ) )
        {
            VSP_error("Level already unselected", DBG_INFO) ;
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
            VSP_error("Using wrong function to select levels", DBG_INFO) ;
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
            VSP_error("Using wrong function to select levels", DBG_INFO) ;
        }
    #endif

    levels->num_bot_levels = n_bot_levels;

    return ;
}

/* Create a deep copy */
VSP_levels* VSP_levels_copy(VSP_levels* levels_in)
{
    VSP_levels *levels_copy = malloc(sizeof(*levels_copy)) ;
    levels_copy->top_levels                = levels_in->top_levels ;
    levels_copy->bot_levels                = levels_in->bot_levels ;
    levels_copy->all_levels                = levels_in->all_levels ;
    levels_copy->no_levels                 = levels_in->no_levels ;
    levels_copy->specify_selected_levels   
        = levels_in->specify_selected_levels ;
    levels_copy->specify_unselected_levels 
        = levels_in->specify_unselected_levels ;

    levels_copy->num_top_levels            = levels_in->num_top_levels ;
    levels_copy->num_bot_levels            = levels_in->num_bot_levels ;

    if(levels_copy->specify_selected_levels)
    {
        levels_copy->levels_selected = VSP_dyn_stk_copy
                                        (
                                            levels_in->levels_selected
                                        ) ;
    }
    else
    {
        levels_copy->levels_selected = NULL ;
    
    }

    if(levels_copy->specify_unselected_levels)
    {
        levels_copy->levels_not_selected = VSP_dyn_stk_copy
                                            (
                                                levels_in->levels_not_selected
                                            ) ;
    }
    else
    {
        levels_copy->levels_not_selected = NULL ;
    }

    return levels_copy ;
}
#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdbool.h>

typedef struct VSP_dyn_stk_struct
{
    /* Public */
    VSPINT *stack     ;
    VSPINT  n_stack   ;

    /* Private */
    VSPINT max_size ;
} VSP_dyn_stk ;


/* Usage:

    First, create a VSP_levels_selection_type variable with the desired
    selection type, and then create a call the constructor of VSP_levels
    and store the result in a variable, called, let's say, levels.

    TOP_LEVELS: Then the top (levels->num_top_levels). Note that if this is 0,
                then no level is selected.

    BOT_LEVELS: Then the bottom (levels->num_bot_levels). Note that if this is
                0, then no level is selected.

    ALL_LEVELS: Then all levels are selected.

    NO_LEVELS: Then no levels are selected.

    SPECIFY_SELECTED_LEVELS: Then the levels_selected stack should be filled
                             with the levels to be selected. This stack is
                             added to using the VSP_levels_select function.

    SPECIFY_UNSELECTED_LEVELS: Then the levels_not_selected stack should be
                               filled with the levels to be unselected. This
                               stack is added to using the VSP_levels_unselect
                               function.
*/
typedef struct VSP_levels_struct
{
    /* Options. Precisely one of the following should be true at any time */
    bool top_levels ;
    bool bot_levels ;
    bool all_levels ;
    bool no_levels ;
    bool specify_selected_levels ;
    bool specify_unselected_levels ;

    /* If 1, only the top (i.e. finest) level. If 0, no level */
    VSPINT num_top_levels ; 
    /* If 1, only the bottom (i.e. coarsest) level. If 0, no level */
    VSPINT num_bot_levels ;

    /* Must be empty unless specify_selected_levels == true */
    VSP_dyn_stk *levels_selected ;

     /* Must be empty unless specify_unselected_levels == true */
    VSP_dyn_stk *levels_not_selected ;

} VSP_levels ;

/*----------------------------------------------------------------------------*/
/* Restart Decision */
/*----------------------------------------------------------------------------*/
/* Choose when to restart solving a level from some new point, like the
   centroid. When constructed, everything is off by default. The user
   must make some selection */




typedef struct VSP_coarsen_dec_struct
{
    VSP_coarsen_dec_type decType ;

    /* Used for useCommunityAtCertainLevels. The levels where comm is ENABLED */
    VSP_levels *levels ;

    /*  If problem->ratioMaxCostNodes > decBoundary, coarsen this problem
        WITHOUT community matching. Otherwise, use community matching.  */
    VSPFLOAT decBounday ;

} VSP_coarsen_dec ;

/*----------------------------------------------------------------------------*/
/* Dynamic Stack */
/*----------------------------------------------------------------------------*/


VSP_dyn_stk *VSP_construct_dyn_stk() ;
void VSP_destruct_dyn_stk(VSP_dyn_stk* ds) ;
void VSP_dyn_stk_push(VSP_dyn_stk* ds, VSPINT to_push) ;
bool VSP_dyn_stk_is_member(VSPINT i, VSP_dyn_stk* ds) ;
VSP_dyn_stk *VSP_dyn_stk_copy(VSP_dyn_stk *stack_in) ;



VSP_coarsen_dec *VSP_construct_coarsen_dec(VSP_coarsen_dec_type decType) ;

#endif
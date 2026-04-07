#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdbool.h>
#include "VSPQP.h"

typedef struct dyn_stk_struct
{
    /* Public */
    VSPINT *stack     ;
    VSPINT  n_stack   ;

    /* Private */
    VSPINT max_size ;
} dyn_stk ;


/* Usage:

    First, create a Levels_selection_type variable with the desired
    selection type, and then create a call the constructor of Levels
    and store the result in a variable, called, let's say, levels.

    TOP_LEVELS: Then the top (levels->num_top_levels). Note that if this is 0,
                then no level is selected.

    BOT_LEVELS: Then the bottom (levels->num_bot_levels). Note that if this is
                0, then no level is selected.

    ALL_LEVELS: Then all levels are selected.

    NO_LEVELS: Then no levels are selected.

    SPECIFY_SELECTED_LEVELS: Then the levels_selected stack should be filled
                             with the levels to be selected. This stack is
                             added to using the Levels_select function.

    SPECIFY_UNSELECTED_LEVELS: Then the levels_not_selected stack should be
                               filled with the levels to be unselected. This
                               stack is added to using the Levels_unselect
                               function.
*/
typedef struct Levels_struct
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
    dyn_stk *levels_selected ;

     /* Must be empty unless specify_unselected_levels == true */
    dyn_stk *levels_not_selected ;

} Levels ;

/*----------------------------------------------------------------------------*/
/* Dynamic Stack */
/*----------------------------------------------------------------------------*/


dyn_stk *construct_dyn_stk() ;
void destruct_dyn_stk(dyn_stk* ds) ;
void dyn_stk_push(dyn_stk* ds, VSPINT to_push) ;
bool dyn_stk_is_member(VSPINT i, dyn_stk* ds) ;
dyn_stk *dyn_stk_copy(dyn_stk *stack_in) ;



#endif
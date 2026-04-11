#ifndef VSPQP_H
#define VSPQP_H

#include "VSPQP_config.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

typedef enum initial_guess_enum
{
    centroid_initial_guess,
    random_initial_guess,
    manual_initial_guess
} initial_guess ;

typedef struct Levels_struct Levels ;

typedef enum
{
    TOP_LEVELS,
    BOT_LEVELS,
    ALL_LEVELS,
    NO_LEVELS,
    SPECIFY_SELECTED_LEVELS,
    SPECIFY_UNSELECTED_LEVELS
} Levels_selection_type ;

Levels* construct_levels
(
    Levels_selection_type selection_type
) ;

typedef struct restart_decision_struct
{
    /* Only one of use_levels, use_rule should be true */
    bool        use_levels ;
    Levels     *levels ; /* Null if we don't use levels */

    /* If the relative change in objective at this level is below rule_param,
       then do a restart */
    bool                   use_rule ;
    VSPINT            restarts_done ;
    VSPINT         restarts_allowed ;
    VSPFLOAT             rule_param ;
    VSPINT     finest_level_allowed ;

} restart_decision ;


typedef struct Options_struct
{
    /* Coarsening 
       If use_coarsening is false, all other coarsening parameters are ignored*/
    bool                    use_coarsening ;
    VSPINT                  max_size_of_coarsest_graph ;

    /* Initial guess */
    initial_guess      initial_guess;
    VSPFLOAT               random_point_weight_proportion ;
    VSPFLOAT               centroidVal ;

    restart_decision *restart_decision ;

    /* Gamma reduction 
       If use_gamma_reduction is false, all other g.r. parameters are ignored */
    bool  use_gamma_reduction ;
    /* Since gamma reduction is much cheaper for small graphs, we may want to
       do it more times for small graphs than for large graphs.               */
    VSPINT num_gr_for_small_graphs ;  /* How many gr to do for small graphs   */
    VSPINT num_gr_for_large_graphs ;  /* How many gr to do for large graphs   */
    VSPINT max_num_vert_small_graph ; /* Max number of vertices in  what is
                                         to be considered a "small" graph     */

    /* Miscellaneous */
    bool     sep_at_each_level ; /* If true, get a separator (i.e. change x, y
                                    to make the penalty term zero) after solving
                                    each level. If false, get a separator only
                                    after the finest level has been solved.   */
    VSPINT      c_is_all_ones  ;
    VSPFLOAT              eps  ; /* Small non-negative number, representing the
                                 tolerance used in various comparisons. As eps
                                 decreases, the solution quality may increase,
                                 but the speed may decrease.                  */


    bool relaxUB ;
    VSPFLOAT ratioUBLimit ;

    /* FM Options */
    VSPINT FM_Limit ;


    /* ====================================================================== */
    // Mountain Climbing Options
    /* ====================================================================== */
    // Tolerance for LP optimality checker
    VSPFLOAT tol_LP_opt_checker ;


    /* ====================================================================== */
    // Simple Coarsening options
    /* ====================================================================== */
    // Options pertaining to the new, simple coarsening

    // Coarsening terminates at a problem P according to the following stopping
    // criteria:
    // 1) If the number of nodes in P is > limit, then 

    // Coarsening terminates at a problem P if
    // 1) P has > limit nodes and P's child is complete, or
    // 2) P has <= limit nodes and the density of P's child is > max_density
    VSPFLOAT max_density     ;  // default = 0.8
    VSPINT   limit_num_nodes ;  // default = 100

    // Enable Hager BM if after HEM and pivot BM, the ratio
    // num coarse nodes / num fine nodes is > compression_threshold
    // To enable Hager BM always, set compression_threshold = 0'
    // To disable Hager BM always, set compression_threshold = 1
    VSPFLOAT compression_threshold ;

} Options;


typedef struct {
    VSPINT  weight_x ;
    VSPINT  weight_y ;
    VSPFLOAT obj_val ;
    VSPINT *x ;
    VSPINT *y ;
} solve_result ;




restart_decision *construct_restart_decision() ;


void Levels_unselect
(
    VSPINT     level_num,
    Levels   *levels
) ;


solve_result *VSP_Solve
(
    VSPINT    n,
    VSPLONG  *Ap,
    VSPINT   *Ai,
    VSPINT   *Ae,      // Temporary argument while I'm refactoring
    Options  *options,
    char     *name     // Optional. Problem name
) ;

void destruct_options(Options* options) ;

void restart_decide_reset(restart_decision *dec) ;

Options* construct_options() ;
void set_default_options(Options* options) ;


long double compute_elapsed_time(struct timespec start, struct timespec end) ;
;

#endif
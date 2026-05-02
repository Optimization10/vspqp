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
    manual_initial_guess,
    bisection_initial_guess
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
    /* Initial guess */
    initial_guess      initial_guess;
    VSPFLOAT               random_point_weight_proportion ;
    VSPFLOAT               centroidVal ;

    restart_decision *restart_decision ;


    /* Since gamma reduction is much cheaper for small graphs, we may want to
       do it more times for small graphs than for large graphs.               */
    VSPINT num_gr_for_small_graphs ;  /* How many gr to do for small graphs   */
    VSPINT num_gr_for_large_graphs ;  /* How many gr to do for large graphs   */
    VSPINT max_num_vert_small_graph ; /* Max number of vertices in  what is
                                         to be considered a "small" graph     */

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
    // Enable Hager BM if after HEM and pivot BM, the ratio
    // num coarse nodes / num fine nodes is > compression_threshold
    // To enable Hager BM always, set compression_threshold = 0'
    // To disable Hager BM always, set compression_threshold = 1
    VSPFLOAT compression_threshold ;
    bool     sort_cols ;
    VSPINT   node_cutoff ;
    VSPFLOAT density_cutoff ;

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
void Print_Options(const Options *opts, const char *path) ;


long double compute_elapsed_time(struct timespec start, struct timespec end) ;
;

#endif
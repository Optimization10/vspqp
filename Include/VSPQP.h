#ifndef VSPQP_H
#define VSPQP_H

#include "VSPQP_config.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>


typedef struct VSP_levels_struct VSP_levels ;
typedef struct VSP_coarsen_dec_struct VSP_coarsen_dec ;


typedef enum VSP_initial_guess_enum
{
    VSP_centroid_initial_guess,
    VSP_random_initial_guess,
    VSP_manual_initial_guess
} VSP_initial_guess ;

typedef struct VSP_restart_decision_struct
{
    /* Only one of use_levels, use_rule should be true */
    bool        use_levels ;
    VSP_levels     *levels ; /* Null if we don't use levels */

    /* If the relative change in objective at this level is below rule_param,
       then do a restart */
    bool                   use_rule ;
    VSPINT            restarts_done ;
    VSPINT         restarts_allowed ;
    VSPFLOAT             rule_param ;
    VSPINT     finest_level_allowed ;

} VSP_restart_decision ;

typedef struct Options_struct
{

    /* Coarsening 
       If use_coarsening is false, all other coarsening parameters are ignored*/
    bool                    use_coarsening ;
    VSPINT                  max_size_of_coarsest_graph ;
    bool                    do_community_match ;


    /* Initial guess */
    VSP_initial_guess      initial_guess;
    VSPFLOAT               random_point_weight_proportion ;
    VSPFLOAT               centroidVal ;



    /* Gamma reduction 
       If use_gamma_reduction is false, all other g.r. parameters are ignored */
    bool  use_gamma_reduction ;
    /* Since gamma reduction is much cheaper for small graphs, we may want to
       do it more times for small graphs than for large graphs.               */
    VSPINT num_gr_for_small_graphs ;  /* How many gr to do for small graphs   */
    VSPINT num_gr_for_large_graphs ;  /* How many gr to do for large graphs   */
    VSPINT max_num_vert_small_graph ; /* Max number of vertices in  what is
                                         to be considered a "small" graph     */

    /* Number of times to do MCA + GR */
    VSPINT num_its_for_small_graphs ;
    VSPINT num_its_for_large_graphs ;

    /* Restart level */
    VSP_restart_decision *restart_decision ;


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
} VSP_solve_result ;


typedef enum
{
    TOP_LEVELS,
    BOT_LEVELS,
    ALL_LEVELS,
    NO_LEVELS,
    SPECIFY_SELECTED_LEVELS,
    SPECIFY_UNSELECTED_LEVELS
} VSP_levels_selection_type ;

VSP_restart_decision *VSP_construct_restart_decision() ;
VSP_levels* VSP_construct_levels
(
    VSP_levels_selection_type selection_type
) ;

void VSP_levels_unselect
(
    VSPINT     level_num,
    VSP_levels   *levels
) ;







VSP_solve_result *VSP_Solve
(
    VSPINT             n,
    VSPLONG          *Ap,
    VSPINT           *Ai,
    VSPINT *          Ae, // Temporary argument while I'm refactoring
    Options *options
) ;

void VSP_destruct_options(Options* options) ;

void VSP_restart_decide_reset(VSP_restart_decision *dec) ;

Options* VSP_construct_options() ;
void VSP_set_default_options(Options* options) ;


long double VSP_compute_elapsed_time(struct timespec start, struct timespec end) ;

#endif
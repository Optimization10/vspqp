#ifndef VSPQP_H
#define VSPQP_H

#include "VSPQP_config.h"
#include <stdbool.h>





typedef struct
{
    /* Public */
    VSPINT *stack     ;
    VSPINT  n_stack   ;

    /* Private */
    VSPINT max_size ;
} VSP_dyn_stk ;

typedef struct
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
typedef struct
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

typedef enum VSP_coarsen_dec_type_enum
{
    CRSN_DEC_useCommAtAllLevels,
    CRSN_DEC_useCommAtNoLevels,
    CRSN_DEC_useCommAtCertainLevels,
    CRSN_DEC_useCommWithBoundary
} VSP_coarsen_dec_type ;

typedef enum VSP_initial_guess_enum
{
    VSP_centroid_initial_guess,
    VSP_random_initial_guess,
    VSP_manual_initial_guess
} VSP_initial_guess ;

typedef struct
{
    VSP_coarsen_dec_type decType ;

    /* Used for useCommunityAtCertainLevels. The levels where comm is ENABLED */
    VSP_levels *levels ;

    /*  If problem->ratioMaxCostNodes > decBoundary, coarsen this problem
        WITHOUT community matching. Otherwise, use community matching.  */
    VSPFLOAT decBounday ;

} VSP_coarsen_dec ;

typedef enum VSP_matching_strategy_enum
{
    HEM    = 0,
    HEMSR  = 1,
    LOW_GAMMA_MATCHING,

} VSP_matching_strategy ;

typedef struct VSP_options_struct
{

    /* Coarsening 
       If use_coarsening is false, all other coarsening parameters are ignored*/
    bool                    use_coarsening ;
    VSPINT                  max_size_of_coarsest_graph ;
    bool                    do_community_match ;
    VSP_matching_strategy   matching_strategy;
    VSP_coarsen_dec        *coarsenDec ;

    // Value by which coarsening should compress the graph
    VSPFLOAT compress_threshold ;

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

} VSP_options;

#endif
#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include<stdbool.h>

typedef struct Stopwatches Stopwatches ;

// #if VSP_INTEGERS_MODE == 0
//     #define VSPINT int ;
//     #define VSPLONG int ;
// #elif VSP_INTEGERS_MODE == 1
//     #define VSPINT int ;
//     #define VSPLONG long long ;
// #elif VSP_INTEGERS_MODE == 2
//     #define VSPINT long long ;
//     #define VSPLONG long long ;
// #else
//     #error "Invalid VSP_INTEGERS_MODE"
// #endif

// typedef double VSPFLOAT ;

typedef enum
{
 VERSION_DEFAULT,
 VERSION_1,
 VERSION_2,
 VERSION_3,
 VERSION_4,
 VERSION_5,
 VERSION_6,
 VERSION_7,
 VERSION_8
} VersionNum ;

typedef enum
{
    STRICT_INCREASE,
    WEAK_INCREASE,
    NO_CHANGE,
    WEAK_DECREASE,
    STRICT_DECREASE,
    UNSPECIFIED_EFFECT
} Obj_Effect ;

typedef struct
{
    /* Public */
    VSPINT *stack     ;
    VSPINT  n_stack   ;

    /* Private */
    VSPINT max_size ;
} VSP_dyn_stk ;

typedef enum
{
    N1_OPTIMAL,            // No feasible, strictly beneficial flip
    N1_WEAKLY_OPTIMAL,     // No feasible, strictly beneficial flip AND there is at least one neutral downflip
    N1_STRICTLY_OPTIMAL,   // No feasible, strictly beneficial flip AND there are zero neutral downflips
    NOT_N1_OPTIMAL         // At least one feasible, strictly beneficial flip
} VSP_N1_Optimality ;

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
typedef enum
{
    TOP_LEVELS,
    BOT_LEVELS,
    ALL_LEVELS,
    NO_LEVELS,
    SPECIFY_SELECTED_LEVELS,
    SPECIFY_UNSELECTED_LEVELS
} VSP_levels_selection_type ;

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

typedef struct {

    VSPINT num_nodes_total ;
    VSPINT num_edges_total ;

    VSPINT num_nodes_x ;
    VSPINT num_nodes_y ;
    VSPINT num_nodes_sep ;

    VSPFLOAT percent_nodes_x ;
    VSPFLOAT percent_nodes_y ;
    VSPFLOAT percent_nodes_sep ;

    VSPFLOAT time ;
} VSP_solution_stats ;


typedef struct {
    VSPINT  weight_x ;
    VSPINT  weight_y ;
    VSPFLOAT obj_val ;
    VSPINT *x ;
    VSPINT *y ;
} VSP_solve_result ;


typedef enum VSP_initial_guess_enum
{
    VSP_centroid_initial_guess,
    VSP_random_initial_guess,
    VSP_manual_initial_guess
} VSP_initial_guess ;


typedef enum VSP_matching_strategy_enum
{
    HEM    = 0,
    HEMSR  = 1,
    LOW_GAMMA_MATCHING,

} VSP_matching_strategy ;



typedef enum VSP_match_type_enum
{
    matchtype_unmatched = 0,
    matchtype_orphan,
    matchtype_standard,
    matchtype_brotherly,
    matchtype_community,

    // Matchtypes for even coarsening
    matchtype_selfMatch,
    matchtype_singleton,
    matchtype_HEM,
    matchtype_3way
} VSP_match_type;



/* Store an almost binary point.                                              */
typedef struct {
    VSPINT    *int_part ; /* Integer part*/
    VSPINT     ind_frac ; /* Index of fractional entry. EMPTY if none.        */
    VSPFLOAT   val_frac ; /* Value of fractional entry. No guarantees
                             if ind_frac is EMPTY.                            */
    VSPINT            n ; /* Number of entries in the point                   */
} VSP_point ;

typedef struct
{
    /*--- Generic Info -------------------------------------------------------*/
    VSPINT n  ;      /* Number of vertices */
    VSPLONG nz ;      /* Twice the number of edges */
    char name[200] ; /* Name of the graph */

    /*--- Resolve from Centroid ----------------------------------------------*/
    /* The levels where we tried re-solving */
    VSP_dyn_stk *resolve_tried ;

    /* Stack containing the levels where re-solving a level
       gave a better point than the current best */
    VSP_dyn_stk *resolve_better ;

    /* Stack containing the initial objective at a level */
    VSP_dyn_stk *objIn  ;

    /* Stack containing the final objective at a level */
    VSP_dyn_stk *objOut ;

    /* Stack containing the value of x'Ay */
    VSP_dyn_stk *numPen ;

    /* Stack containing the sum of all the costs in the separator */
    VSP_dyn_stk *sepIn ;

    /* Stack containing the sum of all the costs in the separator */
    VSP_dyn_stk *sepOut ;

    // The number of flips done in the problem
    long long unsigned num_flips ;

} VSP_stats ;

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

typedef enum
{
    FLIPS_DEC_LEVELS,        /* Specified using a levels struct */
    FLIPS_DEC_NO_COMM,        /* Do flips if the current level was obtained w/o community*/
    FLIPS_DEC_BOUNDARY    /* Do flips if problem->ratioMaxCostNodes > decBoundary */
} VSP_flips_dec_type ;

typedef struct
{
    VSP_flips_dec_type decType ;

    VSP_levels *levels ;
    VSPFLOAT decBoundary ;

} VSP_flips_dec ;

/* Store Az where A is a binary matrix and z has at most one fractional entry
   Entry i in Az is
   int_part[i] + val_frac_z * ( col_frac_z[i] == one ) * (ind_frac_z != EMPTY)*/
typedef struct {
    VSPINT     *int_part ; /* A times the integer part of z                   */
    VSPINT    ind_frac_z ; /* Index of fractional entry in z. EMPTY if none.  */
    VSPFLOAT  val_frac_z ; /* Value of fractional entry in z. No guarantees if
                              if ind_frac_z is EMPTY.                       */
    VSPINT   *col_frac_z ; /* Column of A corresponding to frac index in z. No
                              guarantees if ind_frac_z is EMPTY.              */
    VSPINT           one ;    /* Largest integer in col_frac_v. Represents the
                              1'S in that column. All smaller integers represent
                              zeros.                                          */
    VSPINT             n ; /* Number of entries in Az.                        */
} VSP_Az ;

/* Container for all the information pertaining to x or y */
typedef struct {
    VSP_point         *point ; /* x or y                                      */
    VSP_Az               *Az ; /* Ax or Ay                                    */
    VSPINT            weight ; /* c'x or c'y                                  */
} VSP_state ;


typedef struct
{
    VSPINT  n ;
    VSPLONG nz ;

    VSPINT  *Ai ;
    VSPLONG *Ap ;
    VSPINT  *Ae ;

    VSPINT *c ;
    VSPINT min_vtx_cost ;
    VSPINT max_vtx_cost ;
    VSPINT total_vtx_cost ;
} VSP_Graph ;

typedef enum
{
    UNCOMPRESSED,
    COMPRESSED
} Compression_Status ;

typedef struct VSP_problem_struct
{
    // This is the graph which is being solved
    VSP_Graph             *Graph ;

    char                   *name ;  /* name of the graph                      */

    /* Further problem info */
    VSPINT                     u ;  /* upper balance constraint               */

    /* Info for matching */
    VSPINT             *matching ;
    VSPINT             *matchmap ;
    VSPINT          *invmatchmap ;
    VSP_match_type    *matchtype ;
    VSPINT             singleton ;

    /* Coarser graph size */
    VSPINT                    cn ;
    VSPLONG                   cnz ;

    /* Link to parent problem. NULL if this is finest graph */
    struct VSP_problem_struct* parent_problem ;
    bool was_coarsened_w_comm ;

    /* Link to child problem. NULL if this is the coarsest graph */
    struct VSP_problem_struct* child_problem ;

    /* The number of coarsenings which were done to reach this Graph-> The
       finest problem is at level 0 */
    VSPINT level ;
    VSPINT coarsest_level_num ; /* The level of the coarsest problem */

    /* Re-solve */
    bool resolveTried ;
    bool resolveSuccess ;

    VSPINT ubFinest ;

    /* The ratio of nodes which have maximal cost*/
    VSPFLOAT ratioMaxCostNodes ;

    /* The maximum degree in the graph */
    VSPINT maxDegree ;

    VSP_stats *stats ;

    Stopwatches *Stopwatches ;

} VSP_problem ;

/* The algorithm always performs heavy edge matching and brotherly matching.
   This struct records the user's decision about when to enable community
   matching */

typedef enum VSP_coarsen_dec_type_enum
{
    CRSN_DEC_useCommAtAllLevels,
    CRSN_DEC_useCommAtNoLevels,
    CRSN_DEC_useCommAtCertainLevels,
    CRSN_DEC_useCommWithBoundary
} VSP_coarsen_dec_type ;

typedef struct
{
    VSP_coarsen_dec_type decType ;

    /* Used for useCommunityAtCertainLevels. The levels where comm is ENABLED */
    VSP_levels *levels ;

    /*  If problem->ratioMaxCostNodes > decBoundary, coarsen this problem
        WITHOUT community matching. Otherwise, use community matching.  */
    VSPFLOAT decBounday ;

} VSP_coarsen_dec ;

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

void LoggerInit(bool enable);
void LoggerShutdown(void);
void LoggerStart(VSP_state *Xinit, VSP_state *Yinit, VSP_problem *problem, const char *format, ...);
void LoggerEnd(VSP_state *Xfin, VSP_state *Yfin, VSP_problem *problem, const char *format, ...);

#endif
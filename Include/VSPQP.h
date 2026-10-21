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


typedef enum cover_type_enum
{
    MIN_COST_COVER,
    MIN_CARD_COVER,
    APPROX_MIN_CARD_COVER
} cover_type ;

typedef struct Options_struct
{
    /* Initial guess */
    initial_guess  initial_guess;
    VSPFLOAT       random_point_weight_proportion ;
    VSPFLOAT       centroidVal ;

    restart_decision *restart_decision ;


    /* Since gamma reduction is much cheaper for small graphs, we may want to
       do it more times for small graphs than for large graphs.               */
    VSPINT num_gr_for_small_graphs ;  /* How many gr to do for small graphs   */
    VSPINT num_gr_for_large_graphs ;  /* How many gr to do for large graphs   */
    VSPINT max_num_vert_small_graph ; /* Max number of vertices in  what is
                                         to be considered a "small" graph     */

    /* Gamma reduction only runs on levels at least this coarse. 0 runs it on
       every level. For VSP_Solve_2, branch_shared_levels keeps it in the
       branches and on the shared graph they start from.                      */
    VSPINT gamma_reduction_min_level ;

    bool     tighten_u ;
    VSPFLOAT u_ratio_finest ;
    VSPFLOAT u_ratio_coarsest ;
    bool     ub_based_on_bisection ;

    /* FM Options */
    VSPINT FM_Limit ;


    /* ====================================================================== */
    // Mountain Climbing Options
    /* ====================================================================== */
    // Tolerance for LP optimality checker
    VSPFLOAT tol_LP_opt_checker ;
    int start_mc_both_vars ;


    /* ====================================================================== */
    // Simple Coarsening options
    /* ====================================================================== */
    // Enable Hager BM if after HEM and pivot BM, the ratio
    // num coarse nodes / num fine nodes is > compression_threshold
    // To enable Hager BM always, set compression_threshold = 0'
    // To disable Hager BM always, set compression_threshold = 1
    VSPFLOAT compression_threshold ;
    VSPINT   node_cutoff ;
    VSPFLOAT density_cutoff ;

    // If true, the adjacency lists of every coarse graph are sorted. If false,
    // only a graph handed to vs() (the bisection) is sorted, right before the
    // call, which saves about a third of the coarsening time. Ties in
    // matching, FM and flips follow adjacency order, so results change.
    bool     sort_coarse_graphs ;

    int initial_project ;

    /* ====================================================================== */
    // Vertex Cover Options
    /* ====================================================================== */
    // First, select a cover type:
    // MIN_COST_COVER -- find a min cost cover by solving a max flow problem
    // MIN_CARD_COVER -- find a min cardinality using a max matching
    //
    // APPROX_MIN_CARD_COVER -- use the pricing algorithm to find a cover with
    //                          cardinality <= 2*OPT where OPT is the min
    //                          cardinality among all vertex covers
    cover_type cover_type ;

    // if cover_type = MIN_CARD_COVER, this boolean option controls whether
    //                 we compute a cover starting from each side and select
    //                 the cover with lower cost
    int try_both_covers ;

    /* ====================================================================== */
    // Bisection ngp options
    /* ====================================================================== */

    int bisection_ngp ;

    /* ====================================================================== */
    // FMsubspace options   
    /* ====================================================================== */
    int use_FM_subspace_grad_proj ;
    VSPFLOAT noFM ;

    /* ====================================================================== */
    // Second start options
    /* ====================================================================== */
    // In the second start, if the imbalance is sufficiently large after
    // optimizing at a coarse level, we apply vs() to this point, followed by
    // the optimization routine, and accept this point if it has a (weakly)
    // smaller separator.

    // 0 or 1. Controls if the second start is attempted. Other fields ignored
    // if 0
    int      attempt_second_start ;
    
    // min value of |C(X) - C(Y)| / C(V) required to trigger second start
    VSPFLOAT min_imbal_ratio_second_start ;

    // min level required to trigger second start
    VSPINT   min_level_second_start ;

    /* ====================================================================== */
    // Gradient projection applied to vsp
    /* ====================================================================== */
    // Flag to control if we do vsp_grad_proj
    bool vsp_grad_proj ;
    
    // Number of gradient projection steps applied to vsp
    int  vsp_ngp ; 


    /* ====================================================================== */
    // Hub removal
    /* ====================================================================== */
    double hub_treshold ;

    /* ====================================================================== */
    // Band refinement
    /* ====================================================================== */
    // Optimize runs Band_Refine (a minimum cut on the separator plus 3 BFS
    // layers into each shore) before its first greedy flips, on levels at
    // most band_max_level. EMPTY turns it off. Each call makes at most
    // band_passes improving moves.
    VSPINT band_max_level ;
    VSPINT band_passes ;

    /* ====================================================================== */
    // Branched coarsening (VSP_Solve_2 only)
    /* ====================================================================== */
    // VSP_Solve_2 coarsens the finest graph branch_shared_levels times. From
    // that graph, it builds num_branches independent coarsenings and solves
    // each of them back up to the shared graph. The best of these solutions
    // is refined up to the finest level. Branch 0 matches vertices in index
    // order, like VSP_Solve. The others match them in a random order.
    VSPINT branch_shared_levels ;
    VSPINT num_branches ;

    // VSP_Solve_2 stops coarsening, in the shared levels and in the branches,
    // once a graph has fewer than this many vertices. It is VSP_Solve_2's
    // node_cutoff.
    VSPINT branch_node_cutoff ;

    // Seeds the random visit orders. They use their own generator, so they
    // don't touch the rand() stream that the solvers use.
    VSPINT branch_seed ;

    // How many branches, counting from branch 0, get the full level solve on
    // the levels they solve on their own (their coarsest level up to and
    // including the shared graph). The other branches get a cheap solve there:
    // greedy flips and FM only. If the winner is a cheap branch, it then gets
    // the full level solve on the shared graph. Every finer level always gets
    // the full solve. 0 is the fastest setting.
    VSPINT branch_full_solves ;

} Options ;


typedef struct
{
    VSPINT    weight_x ;
    VSPINT    weight_y ;
    VSPFLOAT  obj_val ;
    VSPINT   *x ;
    VSPINT   *y ;
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
    Options  *options,
    char     *name     // Optional. Problem name
) ;

// Like VSP_Solve, but with METIS-style branched coarsening. See the branch_*
// fields of Options.
solve_result *VSP_Solve_2
(
    VSPINT    n,
    VSPLONG  *Ap,
    VSPINT   *Ai,
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

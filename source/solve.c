/* Thile file contains the VSP_solve function, which is the function the user
   is meant to call. This function handles coarsening the graph, refining the
   graph, and solving the VSP at each level                                   */

#include "Stopwatches.h"
#include "shared_api.h"
#include "solve.h"
#include "FM_api.h"
#include "flips_api.h"
#include "FM_flips_api.h"
#include <time.h>
#include "flips_public_v4.h"
#include "flips_public_v5.h"
#include "gr_flips.h"
#include "logger.h"
#include "globals.h"
#include "coarsen_api.h"
#include "limits.h"



Logger LOG_G ;



void Destroy_Coarse_Problems(VSP_problem *problemFinest)
{
    VSP_problem *cp, *p = problemFinest->child_problem ;
    while(p != NULL)
    {
        cp = p->child_problem ;
        VSP_destruct_problem(p) ;
        p = cp ;
    }
    return ;
}

// Store x, y into result and free remaining memory
void Final_Cleanup
(
    VSP_state *X, VSP_state *Y, VSP_problem *problemFinest, VSP_solve_result *result
)
{
    result->x = X->point->int_part ;
    result->y = Y->point->int_part ;

    VSP_destruct_Az(X->Az) ;
    free(X->point) ;
    free(X) ;

    VSP_destruct_Az(Y->Az) ;
    free(Y->point) ;
    free(Y) ;

    Destroy_Coarse_Problems(problemFinest) ;
}

// Input: an almost binary point. Initial point at this level
// Output: all our improvements algorithms have been applied
void Solve_Level
(
    VSP_state *X, VSP_state *Y, VSP_problem *problem, VSP_options *options
)
{
    printf("\nSolving level: %lld\n", (long long) problem->level) ;
    fflush(stdout) ;
    Optimize(X, Y, problem, options, true, 0) ;
    VSP_relaxed_to_vertex_separator(X, Y, problem) ;
    VSP_ResolveLevelDecideAndDo(X, Y, problem, options) ;
}

VSP_solve_result* VSP_solve
(
    VSP_problem *problem_finest,
    VSP_options        *options
)
{
    /* Preprocessing */
    VSP_check_for_complete_graph(problem_finest) ;
    VSP_preprocess_finest_problem(problem_finest, options) ;
    LoggerInit(false) ;


    /* Coarsening */
    Stopwatches_Start(problem_finest->Stopwatches, "Coarsening") ;
    VSP_problem *problem = VSP_coarsen_all_the_way_down(problem_finest, options) ;
    Stopwatches_Stop(problem_finest->Stopwatches, "Coarsening") ;

    /* Multilevel Solve */
    VSP_state *X, *Y, *Xold, *Yold ;
    for( ; problem != NULL ; problem = problem->parent_problem)
    {
        // printf("\nStarting level %li", CAST(problem->level, long)) ;
        fflush(stdout) ;
        Xold = X ; X = VSP_construct_state(problem->Graph->n) ;
        Yold = Y ; Y = VSP_construct_state(problem->Graph->n) ;

        if( Is_Coarsest_Problem(problem) ) VSP_set_initial_guess(X, Y, problem, options) ;
        else                               Uncoarsen_Point(problem, X, Y, Xold, Yold) ;

        Stopwatches_Start(problem->Stopwatches, "Level") ;
        Solve_Level(X, Y, problem, options) ;
        Stopwatches_Stop(problem->Stopwatches, "Level") ;

        // Print the stopwatch
        Stopwatches_Print(problem->Stopwatches) ;
    }

    /* Cleanup and return */
    VSP_relaxed_to_vertex_separator(X, Y, problem_finest) ;
    VSP_solve_result *result = malloc( sizeof(*result) ) ;
    Final_Cleanup(X, Y, problem_finest, result) ;
    LoggerDestroy() ;

    return result ;
}


/* Determine if we should resolve, and if so, do it */
void VSP_ResolveLevelDecideAndDo
(
    VSP_state             *X,
    VSP_state             *Y,
    VSP_problem     *problem,
    VSP_options     *options
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    bool dec ;
    dec = VSP_restart_decide(options->restart_decision, 0, 0, problem);

    if(dec == false || problem->child_problem == NULL) return ;

    /* Local vars */
    VSPFLOAT obj_best, obj_cent ;

    VSP_state *Xcent = VSP_construct_state(problem->Graph->n) ;
    VSP_state *Ycent = VSP_construct_state(problem->Graph->n) ;

    VSP_set_initial_guess(Xcent, Ycent, problem, options) ;
    Optimize(Xcent, Ycent, problem, options, true, 0) ;
    Force_To_Binary_Point(Xcent, Ycent, problem) ;

    /* If the centroid is better, switch */
    obj_best = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
    obj_cent = VSP_compute_obj(Xcent, Ycent, problem, problem->Graph->max_vtx_cost) ;

    if(obj_cent > obj_best)
    {
        problem->resolveSuccess = true ;

        Swap_Mem(X, Xcent, sizeof( VSP_state )  ) ;
        Swap_Mem(Y, Ycent, sizeof( VSP_state )  ) ;
    }

    VSP_destruct_state(Xcent) ;
    VSP_destruct_state(Ycent) ;
    return ;
}


void VSP_check_for_complete_graph(VSP_problem *problem)
{
    long long n_long  = (long long) problem->Graph->n ;
    long long nz_long = (long long) problem->Graph->nz ;
    long long upper_bnd_nz = n_long * n_long - n_long ;
    if( nz_long == upper_bnd_nz )
    {
        fprintf(stderr, "ERROR: Graph %s is complete\n", problem->name) ;
        exit(1) ;
    }
    else if(nz_long > upper_bnd_nz)
    {
        fprintf(stderr, "ERROR: nz > n^2 - n\n") ;
        exit(1) ;
    }
    return ;
}
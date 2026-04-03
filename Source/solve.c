/* Thile file contains the solve function, which is the function the user
   is meant to call. This function handles coarsening the graph, refining the
   graph, and solving the VSP at each level                                   */

#include "shared_api.h"
#include "solve.h"
#include "FM_api.h"
#include "flips_api.h"
#include <time.h>
#include "flips_public_v4.h"
#include "flips_public_v5.h"
#include "gr_flips.h"
#include "logger.h"
#include "globals.h"
#include "coarsen_api.h"
#include "limits.h"
#include "shared_functions.h"



Logger LOG_G ;



void Destroy_Coarse_Problems(Problem *problemFinest) ;
// {
//     Problem *cp, *p = problemFinest->child_problem ;
//     while(p != NULL)
//     {
//         cp = p->child_problem ;
//         Destroy_Problem(p) ;
//         p = cp ;
//     }
//     return ;
// }

// Store x, y into result and free remaining memory
void Final_Cleanup
(
    State *X, State *Y, Problem *problemFinest, solve_result *result
)
{
    result->x = X->point->int_part ;
    result->y = Y->point->int_part ;

    Destroy_Az(X->Az) ;
    free(X->point) ;
    free(X) ;

    Destroy_Az(Y->Az) ;
    free(Y->point) ;
    free(Y) ;
 
    Destroy_Coarse_Problems(problemFinest) ;
}

// Input: an almost binary point. Initial point at this level
// Output: all our improvements algorithms have been applied
void Solve_Level
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    VSPFLOAT u_ratio = 0.6 ;
    VSPFLOAT s_ratio = 0.47 ;
    VSPFLOAT curr_ratio ;
    VSPINT   L = problem->coarsest_level_num ;

    curr_ratio = (s_ratio - u_ratio) / L + u_ratio ;

    Optimize(X, Y, problem, options, true, 0) ;
    relaxed_to_vertex_separator(X, Y, problem) ;
    ResolveLevelDecideAndDo(X, Y, problem, options) ;
}



/* Determine if we should resolve, and if so, do it */
void ResolveLevelDecideAndDo
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    bool dec ;
    dec = restart_decide(options->restart_decision, 0, 0, problem);

    if(dec == false || problem->child_problem == NULL) return ;

    /* Local vars */
    VSPFLOAT obj_best, obj_cent ;

    State *Xcent = Construct_State(problem->graph->n) ;
    State *Ycent = Construct_State(problem->graph->n) ;

    set_initial_guess(Xcent, Ycent, problem, options) ;
    Optimize(Xcent, Ycent, problem, options, true, 0) ;
    Force_To_Binary_Point(Xcent, Ycent, problem) ;

    /* If the centroid is better, switch */
    obj_best = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
    obj_cent = compute_obj(Xcent, Ycent, problem, problem->graph->max_vtx_cost) ;

    if(obj_cent > obj_best)
    {
        problem->resolveSuccess = true ;

        Swap_Mem(X, Xcent, sizeof( State )  ) ;
        Swap_Mem(Y, Ycent, sizeof( State )  ) ;
    }

    Destroy_State(Xcent) ;
    Destroy_State(Ycent) ;
    return ;
}


void check_for_complete_graph(Problem *problem)
{
    long long n_long  = (long long) problem->graph->n ;
    long long nz_long = (long long) problem->graph->nz ;
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
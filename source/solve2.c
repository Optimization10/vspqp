#include "Stopwatch/Stopwatches.h"
#include "coarsen_api.h"
#include "shared_api.h"
#include "delete_nodes.h"
#include "shared_functions.h"
#include "shared_macros.h"
#include "shared_types.h"
#include <bits/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "simplematch.h"
#include "other_api.h"
#include "solve.h"
#include "mc_api.h"

extern struct timespec  G_START_TIME, G_END_TIME ;
extern long double G_COMPRESS_TIME ;

void SetUpperBound
(
    VSP_problem *problem,
    VSP_options *options
) ;

// Input: an almost binary point. Initial point at this level
// Output: all our improvements algorithms have been applied
void Solve_Level
(
    VSP_state *X, VSP_state *Y, VSP_problem *problem, VSP_options *options
) ;

VSP_solve_result *Solve_Null_Graph
(
    VSP_problem *probFinest
) ;

void Final_Cleanup
(
    VSP_state *X, VSP_state *Y, VSP_problem *problemFinest, VSP_solve_result *result
) ;

void Destroy_Coarse_Problems(VSP_problem *problemFinest) ;

VSP_problem *Problem_Recursive_Coarsen
(
    VSP_problem *probFinest,
    VSP_options    *options
) ;


static void Open_Log_File
(
    long long unsigned count,
    FILE **fp
)
{
    const char *mode = (count == 1) ? "w" : "a";
    *fp = fopen("output/CoarseningLog.txt", mode);
}

static void Log_Problem_Info
(
    FILE *fp,
    VSP_problem *problem
)
{
    fprintf(fp, "Coarsening problem: %s", problem->name) ;
    fflush(fp) ;
}


// Compute total, min, max vertex cost
static void Graph_Compute_Costs
(
    VSP_Graph *Graph
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i, n = Graph->n ;
    VSPINT *c = Graph->c ;
    VSPINT total = 0, max = 0, min = INT_MAX ;

    
    for(i = 0 ; i < n ; i++)
    {
        total += c[i] ;
        if(c[i] > max) max = c[i] ;
        else if(c[i] < min) min = c[i] ;
    }

    Graph->total_vtx_cost = total ;
    Graph->max_vtx_cost = max ;
    Graph->min_vtx_cost = min ;

    return ;
}



VSP_solve_result *Solve_v2
(
    VSP_problem *probFinest,
    VSP_options    *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSP_problem *probCoarsest, *P ;

    /* Preprocessing */
    VSP_check_for_complete_graph(probFinest) ;
    VSP_preprocess_finest_problem(probFinest, options) ;

    #ifdef VSP_DEBUG
        VSP_check_problem(probFinest, options) ;
    #endif

    if(probFinest->Graph->nz == 0)
    {
        VSP_solve_result *ret = Solve_Null_Graph(probFinest) ;
        Destroy_Coarse_Problems(probFinest) ;
        return ret ;
    }



    /*--------------------------------------------------------------------*/
    // Coarsen
    /*--------------------------------------------------------------------*/
    probCoarsest = Problem_Recursive_Coarsen(probFinest, options) ;
    if(probCoarsest == NULL)
    {
        printf("\nExiting for unknown reason. Graph: %s", probFinest->name) ;
        exit(1) ;
    }

    VSP_state *X_l = NULL, *Y_l = NULL, *X_lp1 = NULL , *Y_lp1 = NULL,
              *X_0 = NULL, *Y_0 = NULL ;

    /*------------------------------------------------------------------------*/
    // Solve the coarsest level
    /*------------------------------------------------------------------------*/
    P = probCoarsest ;
    X_lp1 = VSP_construct_state(P->Graph->n) ;
    Y_lp1 = VSP_construct_state(P->Graph->n) ;
    VSP_set_initial_guess(X_lp1, Y_lp1, P, options) ;

    Stopwatches_Start(P->Stopwatches, "Total") ;
    Solve_Level(X_lp1, Y_lp1, P, options) ;
    Stopwatches_Stop(P->Stopwatches, "Total") ;

    // Print the stopwatch
    Stopwatches_Print(P->Stopwatches) ;

    // Invariant: X_l, Y_l holds the solution at level l
    for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
    {
        Stopwatches_Start(P->Stopwatches, "Total") ;

        // Uncoarsen point from level l+1
        Stopwatches_Start(P->Stopwatches, "UncoarsenPoint") ;
        X_l = VSP_construct_state(P->Graph->n) ;
        Y_l = VSP_construct_state(P->Graph->n) ;
        Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;
        Stopwatches_Stop(P->Stopwatches, "UncoarsenPoint") ;

        // Solve level l
        Solve_Level(X_l, Y_l, P, options) ;
        VSP_relaxed_to_vertex_separator(X_l, Y_l, P) ;

        X_lp1 = X_l ;
        Y_lp1 = Y_l ;

        Stopwatches_Stop(P->Stopwatches, "Total") ;

        /* ================================================================== */
        // Check objective
        /* ================================================================== */
        if(X_l->point->ind_frac != EMPTY || Y_l->point->ind_frac != EMPTY)
        {
            printf("\nThere is a fractional entry") ;
            fflush(stdout) ;
        }
        VSPINT cost_shores = 0 ;
        VSPINT cost_graph  = 0 ;
        VSPINT xAy         = 0 ;
        VSPINT *x          = X_l->point->int_part ;
        VSPINT *y          = Y_l->point->int_part ;
        VSPINT *Ay         = Y_l->Az->int_part ;
        VSPINT *c          = P->Graph->c ;
        for(VSPINT i = 0 ; i < P->Graph->n ; i++)
        {
            cost_shores += c[i] * x[i] ;
            cost_shores += c[i] * y[i] ;
            cost_graph  += c[i] ;
            xAy         += x[i] * Ay[i] ;
        }

        if(xAy == 0)
        {
            printf("\nLevel: %i  Sep Cost: %i", P->level, cost_graph - cost_shores) ;
            fflush(stdout) ;
        }
        else
        {
            printf("\nPenalty term not zero") ;
            fflush(stdout) ;
        }


        /* ================================================================== */

        // Print the stopwatch and number of flips
        Stopwatches_Print(P->Stopwatches) ;
        Print_Num_Flips(P) ;
    }

    // X_0, Y_0 is the solution at the finest, uncompressed level
    X_0 = X_lp1 ;
    Y_0 = Y_lp1 ;
    

    // Cleanup
    VSP_relaxed_to_vertex_separator(X_0, Y_0, probFinest) ;
    VSP_solve_result *result = malloc( sizeof(*result) ) ;
    Final_Cleanup(X_0, Y_0, probFinest, result) ;

    return result ;

}







// Solve the VSP in the case that the graph has no edges
// Input: The graph has no edges, all vertex costs are 1
// Output: 
VSP_solve_result *Solve_Null_Graph
(
    VSP_problem *probFinest
)
{
    VSPINT i,
    n    = probFinest->Graph->n,
    u    = probFinest->u ;

    VSPINT M    = MIN(u, n / 2) ;


    VSP_solve_result *result = malloc( sizeof(*result) ) ;
    result->x = calloc(n, sizeof(VSPINT)) ;
    result->y = calloc(n, sizeof(VSPINT)) ;

    if( u < n/2 )
    {
        for(i = 0 ; i < u ; i++)
        {
            result->x[  i  ] = 1 ;
            result->y[u + i] = 1 ;
        }
        result->weight_x = u ;
        result->weight_y = u ;
    }
    else // u >= n/2
    {
        for(i = 0 ; i < n/2 ; i++)
        {
            result->x[  i  ] = 1 ;
            result->y[n/2 + i] = 1 ;
            result->weight_x = n/2 ;
            result->weight_y = n/2 ;
        }

        if(n % 2)
        {
            result->y[n-1] = 1 ;
            result->weight_y += 1 ;
        }
    }
    result->obj_val  = result->weight_x + result->weight_y ;
    Destroy_Coarse_Problems(probFinest) ;
    return result ;
}
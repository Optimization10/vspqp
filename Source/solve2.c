#include "VSPQP.h"
#include "core.h"
#include "coarsen_api.h"
#include "shared_api.h"
#include "delete_nodes.h"
#include "shared_functions.h"

#include <bits/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "other_api.h"
#include "solve.h"
#include "mc_api.h"
#include "coarsen_api.h"

extern struct timespec  G_START_TIME, G_END_TIME ;
extern long double G_COMPRESS_TIME ;

const Problem *G_P = NULL ;

void SetUpperBound
(
    Problem *problem,
    Options *options
) ;

static void Set_Global_Problem(Problem *P)
{
    G_P = P ;
    return ;
}

// Input: an almost binary point. Initial point at this level
// Output: all our improvements algorithms have been applied
void Solve_Level
(
    State *X, State *Y, Problem *problem, Options *options
) ;

solve_result *Solve_Null_Graph
(
    Problem *probFinest
) ;

void Final_Cleanup
(
    State *X, State *Y, Problem *problemFinest, solve_result *result
) ;

void Destroy_Coarse_Problems(Problem *problemFinest)
{
    Problem *cp, *p = problemFinest->child_problem ;
    while(p != NULL)
    {
        cp = p->child_problem ;
        destruct_problem(p) ;
        p = cp ;
    }
    return ;
}


void set_default_problem_data(Problem* problem) ;

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
    Problem *problem
)
{
    fprintf(fp, "Coarsening problem: %s", problem->name) ;
    fflush(fp) ;
}


// Compute total, min, max vertex cost
static void Graph_Compute_Costs
(
    Graph *Graph
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



solve_result *Solve
(
    VSPINT             n,
    VSPLONG          *Ap,
    VSPINT           *Ai,
    VSPINT *          Ae, // Temporary argument while I'm refactoring
    Options *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    Problem *probFinest, *probCoarsest, *P ;

    probFinest = construct_problem(n, Ap[n], Ai, Ap, Ae, options) ;
    set_default_problem_data(probFinest) ;

    /* Preprocessing */
    construct_problem(n, Ap[n], Ai, Ap, NULL, options) ;
    check_for_complete_graph(probFinest) ;
    preprocess_finest_problem(probFinest, options) ;

    #ifdef DEBUG
        check_problem(probFinest, options) ;
    #endif

    if(probFinest->Graph->nz == 0)
    {
        solve_result *ret = Solve_Null_Graph(probFinest) ;
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

    State *X_l = NULL, *Y_l = NULL, *X_lp1 = NULL , *Y_lp1 = NULL,
              *X_0 = NULL, *Y_0 = NULL ;

    /*------------------------------------------------------------------------*/
    // Solve the coarsest level
    /*------------------------------------------------------------------------*/
    P = probCoarsest ;
    X_lp1 = construct_state(P->Graph->n) ;
    Y_lp1 = construct_state(P->Graph->n) ;
    set_initial_guess(X_lp1, Y_lp1, P, options) ;

    Set_Global_Problem(P) ;

    Solve_Level(X_lp1, Y_lp1, P, options) ;


    // Invariant: X_l, Y_l holds the solution at level l
    for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
    {
        // printf("u: %f\n", (double) P->u / P->Graph->total_vtx_cost) ;
        // fflush(stdout) ;

        Set_Global_Problem(P) ;



        // Uncoarsen point from level l+1

        X_l = construct_state(P->Graph->n) ;
        Y_l = construct_state(P->Graph->n) ;
        Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;


        // Solve level l
        Solve_Level(X_l, Y_l, P, options) ;
        relaxed_to_vertex_separator(X_l, Y_l, P) ;

        X_lp1 = X_l ;
        Y_lp1 = Y_l ;
    }

    // X_0, Y_0 is the solution at the finest, uncompressed level
    X_0 = X_lp1 ;
    Y_0 = Y_lp1 ;
    

    // Cleanup
    relaxed_to_vertex_separator(X_0, Y_0, probFinest) ;
    solve_result *result = malloc( sizeof(*result) ) ;
    Final_Cleanup(X_0, Y_0, probFinest, result) ;

    Set_Global_Problem(NULL) ;
    return result ;
}







// Solve the VSP in the case that the graph has no edges
// Input: The graph has no edges, all vertex costs are 1
// Output: 
solve_result *Solve_Null_Graph
(
    Problem *probFinest
)
{
    VSPINT i,
    n    = probFinest->Graph->n,
    u    = probFinest->u ;

    VSPINT M    = MIN(u, n / 2) ;


    solve_result *result = malloc( sizeof(*result) ) ;
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
#include "Stopwatch/Stopwatches.h"
#include "coarsen_api.h"
#include "shared_api.h"
#include "delete_nodes.h"
#include "shared_functions.h"
#include "shared_macros.h"
#include "shared_types.h"
#include <bits/time.h>
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

static VSP_problem *Problem_Recursive_Coarsen
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


// static void Log_Uncompressed_Graph
// (
//     FILE *fp,
//     VSP_problem *problem
// )
// {
//     fprintf(fp, "\n\n    Starting Level:  %lli", CAST(problem->level, long long));
//     fprintf
//     (
//         fp, "\n    Uncompressed graph:  n = %lli  nz = %lli  density = %f",
//         CAST(problem->Graph->n, long long),
//         CAST(problem->Graph->nz, long long),
//         Graph_Density(problem->Graph->n, problem->Graph->nz)
//     ) ;
//     fflush(fp) ;
// }

// static void Log_Compressed_Graph
// (
//     FILE *fp,
//     VSP_problem *problem
// )
// {
//     fprintf
//     (
//         fp, "\n    Compressed graph:  n = %lli  nz = %lli  density = %f",
//         CAST(problem->Graph->n , long long),
//         CAST(problem->Graph->nz, long long),
//         Graph_Density(problem->Graph->n, problem->Graph->nz)
//     ) ;

//     VSPINT perfectCost = pow_ints(2, problem->level) ;
//     double numPerfectCost = 0 ;
//     for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
//     {
//         if(problem->Graph->c[i] == perfectCost) numPerfectCost++ ;
//     }
//     if(numPerfectCost)
//     {
//         fprintf(fp, "\n    Ratio perfect vertices in compressed graph: %f", numPerfectCost / problem->Graph->n) ;
//         fflush(fp) ;
//         return ;
//     }

//     double numMaxCost = 0 ;
//     for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
//     {
//         if(problem->Graph->c[i] == problem->Graph->max_vtx_cost) numMaxCost++ ;
//     }
//     fprintf
//     (
//         fp, "\n    Zero perfect vertices in compressed Graph-> Highest cost: %lli. \n        Ratio of verts with highest cost: %f",
//         CAST(problem->Graph->max_vtx_cost, long long),
//         numMaxCost / problem->Graph->n
//     ) ;
//     fflush(fp) ;

// }

// static void Log_Coarsened_Graph
// (
//     FILE *fp,
//     VSP_problem *problem
// )
// {
//     fprintf
//     (
//         fp, "\n    Coarsened Graph:  n = %lli  nz = %lli  density = %f",
//         CAST(problem->Graph->n, long long),
//         CAST(problem->Graph->nz, long long),
//         Graph_Density(problem->Graph->n, problem->Graph->nz)
//     );
//     fprintf(fp, "\n    Coarsening ratio: %f", problem->Graph->n / CAST(problem->parent_problem->Graph->n, double) ) ;
//     fflush(fp) ;
// }

// static void Set_Compression_Status(VSP_problem *P, Compression_Status newStatus)
// {
//     // Compression_Status currentStatus =  ;
//     if(newStatus == COMPRESSED)
//     {
//         P->Graph = P->Graph ;
//         P->compression_status = COMPRESSED ;
//     }
//     else
//     {
//         P->Graph = P->Graph ;
//         P->compression_status = UNCOMPRESSED ;   
//     }
//     return ;
// }

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

// VSP_solve_result *Solve_v2
// (
//     VSP_problem *probFinest,
//     VSP_options    *options
// )
// {
//     static long long unsigned count = 0 ;
//     count++ ;

//     VSP_problem *probCoarsest, *P ;

//     /* Preprocessing */
//     VSP_check_for_complete_graph(probFinest) ;
//     VSP_preprocess_finest_problem(probFinest, options) ;

//     #ifdef VSP_DEBUG
//         VSP_check_problem(probFinest, options) ;
//     #endif

//     if(probFinest->Graph->nz == 0)
//     {
//         VSP_solve_result *ret = Solve_Null_Graph(probFinest) ;
//         Destroy_Coarse_Problems(probFinest) ;
//         return ret ;
//     }



//     /*--------------------------------------------------------------------*/
//     // Coarsen
//     /*--------------------------------------------------------------------*/
//     probCoarsest = Recursive_Coarsen_And_Compress(probFinest, options) ;
//     if(probCoarsest == NULL)
//     {
//         printf("\nExiting for unknown reason. Graph: %s", probFinest->name) ;
//         exit(1) ;
//     }

//     // Set total vertex cost
//     P = probCoarsest ;
//     while(P != NULL)
//     {
//         Graph_Compute_Costs( &(P->Graph) ) ;
//         Graph_Compute_Costs( &(P->Graph) ) ;
//         Graph_Compute_Costs( &(P->Graph) ) ;
//         P = P->parent_problem ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Solve
//     /*------------------------------------------------------------------------*/
//     // P    = problem at level l
//     // (Xlp1, Ylp1) = final solution at level l+1 (solution to uncomp problem)
//     // (Xl, Yl)     = initial solution at level l
//     // 
//     VSP_state *X_l, *Y_l, *X_lp1, *Y_lp1 ;


//     /*------------------------------------------------------------------------*/
//     // Solve the coarsest level
//     /*------------------------------------------------------------------------*/
//     P = probCoarsest ;

//     // Solve compressed version
//     Set_Compression_Status(P, COMPRESSED) ;

//     X_l = VSP_construct_state(P->Graph->n) ;
//     Y_l = VSP_construct_state(P->Graph->n) ;
//     VSP_set_initial_guess(X_l, Y_l, P, options) ;
    
    
//     Solve_Level(X_l, Y_l, P, options) ;

//     // Project to uncompressed problem
//     clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;
//     Set_Compression_Status(P, UNCOMPRESSED) ;
//     X_lp1 = Uncompress_State
//     (
//         &(P->Graph),
//         P->newIndex,
//         X_l,
//         P->Graph->Ai,
//         P->Graph->Ap,
//         P->Graph->n
//     ) ;

//     Y_lp1 = Uncompress_State
//     (
//         &(P->Graph),
//         P->newIndex,
//         Y_l,
//         P->Graph->Ai,
//         P->Graph->Ap,
//         P->Graph->n
//     ) ;

//     clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
//     G_COMPRESS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;

//     VSP_destruct_state(X_l) ;
//     VSP_destruct_state(Y_l) ;

//     VSPINT loop = 0 ;
//     for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
//     {
//         loop++ ;
//         // Uncoarsen point from level l+1
//         Set_Compression_Status(P, COMPRESSED) ;
//         Set_Compression_Status(P->child_problem, UNCOMPRESSED) ;
//         X_l = VSP_construct_state(P->Graph->n) ;
//         Y_l = VSP_construct_state(P->Graph->n) ;
//         Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;
//         X_lp1 = NULL ;
//         Y_lp1 = NULL ;

//         // Solve the compressed problem
//         Solve_Level(X_l, Y_l, P, options) ;
        

//         // Uncompress the solution
//         clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;
//         Set_Compression_Status(P, UNCOMPRESSED) ;
//         X_lp1 = Uncompress_State
//         (
//             &(P->Graph),
//             P->newIndex,
//             X_l,
//             P->Graph->Ai,
//             P->Graph->Ap,
//             P->Graph->n
//         ) ;

//         Y_lp1 = Uncompress_State
//         (
//             &(P->Graph),
//             P->newIndex,
//             Y_l,
//             P->Graph->Ai,
//             P->Graph->Ap,
//             P->Graph->n
//         ) ;

//         VSP_relaxed_to_vertex_separator(X_lp1, Y_lp1, P) ;

//         clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
//         G_COMPRESS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;


//         VSP_destruct_state(X_l) ;
//         VSP_destruct_state(Y_l) ;
//         X_l = NULL ;
//         Y_l = NULL ;

//     }

//     // X_0, Y_0 is the solution at the finest, uncompressed level
//     VSP_state *X_0 = X_lp1 ;
//     VSP_state *Y_0 = Y_lp1 ;
    

//     // Cleanup
//     VSP_relaxed_to_vertex_separator(X_0, Y_0, probFinest) ;
//     VSP_solve_result *result = malloc( sizeof(*result) ) ;
//     Final_Cleanup(X_0, Y_0, probFinest, result) ;

//     return result ;
// }




// problem = a problem at a certain level
// nOrg = num vertices in problem P
// For vertex i of Porg, keep[i] = 1 if i is not a singleton, 0 if it is
// For vertex i of Porg, newIndex[i] = node in Pcmp corresponding to i if
// keep[i] = 1, and EMPTY otherwise
// VSP_solve_result *Solve_v2
// (
//     VSP_problem *probFinest,
//     VSP_options    *options
// )
// {
//     static long long unsigned count = 0 ;
//     count++ ;

//     VSP_problem *probCoarsest, *P ;

//     /* Preprocessing */
//     VSP_check_for_complete_graph(probFinest) ;
//     VSP_preprocess_finest_problem(probFinest, options) ;

//     #ifdef VSP_DEBUG
//         VSP_check_problem(probFinest, options) ;
//     #endif

//     if(probFinest->Graph->nz == 0)
//     {
//         VSP_solve_result *ret = Solve_Null_Graph(probFinest) ;
//         Destroy_Coarse_Problems(probFinest) ;
//         return ret ;
//     }



//     /*--------------------------------------------------------------------*/
//     // Coarsen
//     /*--------------------------------------------------------------------*/
//     probCoarsest = Recursive_Coarsen_And_Compress(probFinest, options) ;
//     if(probCoarsest == NULL)
//     {
//         printf("\nExiting for unknown reason. Graph: %s", probFinest->name) ;
//         exit(1) ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Solve
//     /*------------------------------------------------------------------------*/
//     // P    = problem at level l
//     // (Xlp1, Ylp1) = final solution at level l+1 (solution to uncomp problem)
//     // (Xl, Yl)     = initial solution at level l
//     // 
//     VSP_state *X_l, *Y_l, *X_lp1, *Y_lp1 ;


//     /*------------------------------------------------------------------------*/
//     // Solve the coarsest level
//     /*------------------------------------------------------------------------*/
//     P = probCoarsest ;

//     // Solve compressed version
//     Set_Compression_Status(P, COMPRESSED) ;

//     X_l = VSP_construct_state(P->Graph->n) ;
//     Y_l = VSP_construct_state(P->Graph->n) ;
//     VSP_set_initial_guess(X_l, Y_l, P, options) ;
    
//     Solve_Level(X_l, Y_l, P, options) ;

//     // Project to uncompressed problem
//     clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;
//     Set_Compression_Status(P, UNCOMPRESSED) ;
//     X_lp1 = Uncompress_State
//     (
//         &(P->Graph),
//         P->newIndex,
//         X_l,
//         P->Graph->Ai,
//         P->Graph->Ap,
//         P->Graph->n
//     ) ;

//     Y_lp1 = Uncompress_State
//     (
//         &(P->Graph),
//         P->newIndex,
//         Y_l,
//         P->Graph->Ai,
//         P->Graph->Ap,
//         P->Graph->n
//     ) ;

//     clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
//     G_COMPRESS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;

//     VSP_destruct_state(X_l) ;
//     VSP_destruct_state(Y_l) ;

//     for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
//     {
//         // Uncoarsen point from level l+1
//         Set_Compression_Status(P, COMPRESSED) ;
//         Set_Compression_Status(P->child_problem, UNCOMPRESSED) ;
//         X_l = VSP_construct_state(P->Graph->n) ;
//         Y_l = VSP_construct_state(P->Graph->n) ;
//         Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;
//         X_lp1 = NULL ;
//         Y_lp1 = NULL ;

//         // Solve the compressed problem
//         Solve_Level(X_l, Y_l, P, options) ;

//         // Uncompress the solution
//         clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;
//         Set_Compression_Status(P, UNCOMPRESSED) ;
//         X_lp1 = Uncompress_State
//         (
//             &(P->Graph),
//             P->newIndex,
//             X_l,
//             P->Graph->Ai,
//             P->Graph->Ap,
//             P->Graph->n
//         ) ;

//         Y_lp1 = Uncompress_State
//         (
//             &(P->Graph),
//             P->newIndex,
//             Y_l,
//             P->Graph->Ai,
//             P->Graph->Ap,
//             P->Graph->n
//         ) ;

//         VSP_relaxed_to_vertex_separator(X_lp1, Y_lp1, P) ;

//         clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
//         G_COMPRESS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;


//         VSP_destruct_state(X_l) ;
//         VSP_destruct_state(Y_l) ;
//         X_l = NULL ;
//         Y_l = NULL ;

//     }

//     // X_0, Y_0 is the solution at the finest, uncompressed level
//     VSP_state *X_0 = X_lp1 ;
//     VSP_state *Y_0 = Y_lp1 ;
    

//     // Cleanup
//     VSP_relaxed_to_vertex_separator(X_0, Y_0, probFinest) ;
//     VSP_solve_result *result = malloc( sizeof(*result) ) ;
//     Final_Cleanup(X_0, Y_0, probFinest, result) ;

//     return result ;
// }

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
    Solve_Level(X_lp1, Y_lp1, P, options) ;

    // Print the stopwatch
    Stopwatches_Print(P->Stopwatches) ;

    // Invariant: X_l, Y_l holds the solution at level l
    for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
    {
        // Uncoarsen point from level l+1
        X_l = VSP_construct_state(P->Graph->n) ;
        Y_l = VSP_construct_state(P->Graph->n) ;
        Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;

        // Solve level l
        Solve_Level(X_l, Y_l, P, options) ;
        VSP_relaxed_to_vertex_separator(X_l, Y_l, P) ;

        X_lp1 = X_l ;
        Y_lp1 = Y_l ;

        // Print the stopwatch
        Stopwatches_Print(P->Stopwatches) ;
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

/* =============================================================================
Compress, Coarsen, Repeat

Input Assumptions:
- Finest level is NOT all singletons

Output Guarantees:
The returned problem is the coarsest level such that:
- Max vtx cost is <= u (why? So that if the only positive entry in z is
  fractional, z can be projected to a feasible binary point by rounding the
  fractional component up to 1)
- 


============================================================================= */
// static VSP_problem *Problem_Recursive_Coarsen
// (
//     VSP_problem *probFinest,
//     VSP_options    *options
// )
// {
//     static long long unsigned count = 0 ;
//     count++ ;

//     FILE *fp;
//     Open_Log_File(count, &fp);

//     VSP_problem *problem, *problemCoarse ;
//     VSPINT i, nFinest, nOrg, *AiOrg, *keep ;
//     VSPLONG *ApOrg ;
//     VSP_Graph *GraphOrg, GraphTemp, *GraphCmp ;
    
//     nFinest = probFinest->Graph->n ;

//     bool stoppingCriterion ;

//     problem = probFinest ;
//     Log_Problem_Info(fp, problem);
//     while( true )
//     {
//         nOrg  = problem->Graph->n  ;
//         ApOrg = problem->Graph->Ap ;
        
//         Log_Uncompressed_Graph(fp, problem) ;
//         /*--------------------------------------------------------------------*/
//         // Compress problem, removing singletons
//         /*--------------------------------------------------------------------*/
//         // Determine the vertices to keep
//         keep = problem->keep ;
//         VSPINT nkeep = 0 ;
//         for(i = 0 ; i < nOrg ; i++)
//         {
//             if(ApOrg[i+1] - ApOrg[i])
//             {
//                 keep[i] = 1 ;
//                 nkeep++ ;
//             }
//             else keep[i] = 0 ;
//         }


//         G_COMPRESS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;

//         bool all_singletons = (GraphCmp == NULL) ;
//         if(GraphCmp != NULL)
//         {
//             problem->Graph = *GraphCmp ;
//             free(GraphCmp) ;
//             GraphCmp = NULL ;
//             Log_Compressed_Graph(fp, problem) ;
//             fflush(stdout) ;
//         }

//         /*--------------------------------------------------------------------*/
//         // Stopping Criteria
//         /*--------------------------------------------------------------------*/
//         stoppingCriterion = false ;

//         // If the graph is made entirely of singletons and compression
//         // returned NULL
//         if(all_singletons)
//         {
//             // Go to the previous graph
//             problem = problem->parent_problem ;
//             VSP_destruct_problem(problem->child_problem) ;
//             stoppingCriterion = true ;   
//         }

//         if(problem->Graph->total_vtx_cost < problem->u)
//         {
//             stoppingCriterion = true ;
//         }

//         if(problem->Graph->n < 10) stoppingCriterion = true ;

//         if(Is_Complete_Graph(problem->Graph->n, problem->Graph->nz))
//         {
//             if(problem->level == 0) return probFinest ; 
            
//             // Go to the previous graph
//             problem = problem->parent_problem ;
//             VSP_destruct_problem(problem->child_problem) ;
//             stoppingCriterion = true ;
//         }

//         if(problem->Graph->max_vtx_cost > problem->u)
//         {
//             // Go to the previous graph
//             problem = problem->parent_problem ;
//             VSP_destruct_problem(problem->child_problem) ;
//             stoppingCriterion = true ;   
//         }

//         if(stoppingCriterion) break ;

//         /*--------------------------------------------------------------------*/
//         // Coarsen the compressed problem
//         /*--------------------------------------------------------------------*/
//         // Match the compressed problem
//         Problem_Simple_Match
//         (
//             &(problem->Graph),
//             problem->level,
//             problem->matching,
//             problem->matchmap,
//             problem->invmatchmap,
//             problem->matchtype,
//             &(problem->cn)
//         ) ;

//         // Coarsen the compressed problem
//         problemCoarse = Problem_Coarsen
//         (
//             &(problem->Graph),
//             problem->cn,
//             problem->matchmap,
//             problem->invmatchmap,
//             problem->matching,
//             options,
//             problem->level,
//             problem->name,
//             problem->l,
//             problem->u,
//             problem->ubFinest
//         ) ;
//         problemCoarse->level = problem->level + 1 ;

//         // Set the total graph weight for the uncompressed problem
//         problemCoarse->Graph->total_vtx_cost = problem->Graph->total_vtx_cost ;

//         problemCoarse->Graph = problemCoarse->Graph ;

//         // Link the child and parent problems
//         problemCoarse->parent_problem = problem ;
//         problem      ->child_problem  = problemCoarse ;
//         problem = problemCoarse ;       
//         problem->Graph = problem->Graph ;

//         Log_Coarsened_Graph(fp, problem) ;        
//     }

//     fprintf(fp, "\n\n") ;
//     fclose(fp) ;

//     // In case we went back a graph
//     problem->child_problem = NULL ;

//     // Set coarsest level num
//     VSPINT coarsest_level_num = problem->level ;

//     VSP_problem *P = problem ;
//     while(P != NULL)
//     {
//         P->coarsest_level_num = coarsest_level_num ;
//         P = P->parent_problem ;
//     }
    
//     /* Set upper bounds */
//     // TODO: Doing things incorrectly here. When lowering the upper balance
//     // constraint at coarser levels, we risk having a vertex whose cost is
//     // greater than the upper balance constraint. This may lead to not being
//     // able to round an almost binary point with a fractional component to a
//     // fully binary point
//     VSP_problem *p = problem ;
//     while(p != NULL)
//     {
//        SetUpperBound(p, options) ;
//     //    printf("Upper bound: %i      Level: %i\n", p->u, p->level) ;
//        p = p->parent_problem ;
//     }

//     return problem ;
// }





VSP_problem *Problem_Recursive_Coarsen
(
    VSP_problem *probFinest,
    VSP_options    *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT n_finest = probFinest->Graph->n ;
    VSP_problem *P, *coarse_P, *coarsest_P ;

    P = probFinest ;
    bool stop = false ;

    // After the loop, P will be the coarsest problem
    VSPINT l = 0 ; // Invariant: l = # of coarsenings
    while( 1 )
    {
        // Match P
        Graph_Simple_Match
        (
            P->Graph,
            P->level,
            P->matching,
            P->matchmap,
            P->invmatchmap,
            P->matchtype,
            &(P->cn)
        ) ;

        // Coarsen P
        coarse_P = Problem_Coarsen(P, options) ;

        // Set the level
        coarse_P->level = P->level + 1 ;

        // Set the total graph weight
        coarse_P->Graph->total_vtx_cost = P->Graph->total_vtx_cost ;

        // Link the child and parent problems
        coarse_P->parent_problem = P ;
        P       ->child_problem  = coarse_P ;    

        // Decide stopping criteria
        if(coarse_P->Graph->n < 10)
        {
            coarsest_P = coarse_P ;
            break ;
        }

        if(Is_Complete_Graph(coarse_P->Graph->n, coarse_P->Graph->nz))
        {
            // Destroy the problem corresponding to a complete graph
            VSP_destruct_problem(coarse_P) ;
            coarsest_P = P ;
            break ;
        }
        
        P = coarse_P ;
        l++ ;
    }

    coarsest_P->child_problem = NULL ;
    return coarsest_P ; // P is coarse problem now and is being returned ####################### Since we hit second stopping criterion
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
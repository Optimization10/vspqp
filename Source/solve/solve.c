#include "solve_internal.h"
#include "VSPQP.h"
#include "core.h"
#include "coarsen_api.h"
#include "utils.h"
#include <string.h>

extern struct timespec  G_START_TIME, G_END_TIME ;
extern long double G_COMPRESS_TIME ;


bool restart_decide
(
    restart_decision           *dec,
    VSPFLOAT                  objIn, /* If not using rule, this is ignored */
    VSPFLOAT                 objOut, /* If not using rule, this is ignored */
    Problem                *problem  /* If using levels, this is ignored   */
) ;



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

    Optimize(X, Y, problem, options) ;
    relaxed_to_vertex_separator(X, Y, problem) ;
    ResolveLevelDecideAndDo(X, Y, problem, options) ;
}

void Destroy_Coarse_Problems(Problem *problemFinest)
{
    Problem *cp, *p = problemFinest->child_problem ;
    while(p != NULL)
    {
        cp = p->child_problem ;
        Destroy_Problem(p) ;
        p = cp ;
    }
    return ;
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
    n    = probFinest->graph->n,
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

solve_result *Package_Result(State *X, State *Y, Problem *P)
{
    VSPINT n = P->graph->n;
    solve_result *result = malloc(sizeof(*result));

    result->x = malloc(n * sizeof(VSPINT));
    result->y = malloc(n * sizeof(VSPINT));
    memcpy(result->x, X->point->int_part, n * sizeof(VSPINT));
    memcpy(result->y, Y->point->int_part, n * sizeof(VSPINT));

    result->weight_x = X->weight;
    result->weight_y = Y->weight;
    result->obj_val  = result->weight_x + result->weight_y;

    return result;
}



/* Set the cost of each vector to 1, each edge weight to 1, and
   the balance constraints as 60% of the total number of vertices             */
void set_default_problem_data(Problem* problem) {
    VSPLONG i ;
    VSPINT n = problem->graph->n ;
    VSPLONG nz = problem->graph->nz ;

    /* Set cost of each vertex to 1 */
    VSPINT* c = problem->graph->c ;
    for(i = 0 ; i < n ; i++) c[i] = 1 ;

    /* Set the edge weights to 1 */
    VSPINT* Ae = problem->graph->Ae ;
    for(i = 0 ; i < nz ; i++) Ae[i] = 1 ;

    /* Set the upper balance constraint to 60%, and lower constraint to 1 */
    problem->u = (VSPINT) (problem->graph->n * 0.6);

    problem->ubFinest = (VSPINT) (problem->graph->n * 0.6) ;

    return ;
}

/* Compute alpha, gamma, sum of all costs, and set the is_original flag to
   true */
void    preprocess_finest_problem
(
    Problem* problem,
    Options* options
)
{
    VSPINT node, num_nodes, c_node, alpha, gamma, total_graph_weight ;
    VSPINT* c ;
    num_nodes = problem->graph->n ;

    /* Finest problem has no parent. Lets us know when to stop refining */
    problem->parent_problem = NULL ;
    problem->level = 0 ;

    /* Compute the maximum degree */
    VSPINT maxDeg = 0 ;
    VSPLONG *Ap    = problem->graph->Ap ;
    for(node = 0 ; node < num_nodes ; node++)
    {
        if( Ap[node + 1] - Ap[node] > maxDeg )
        {
            maxDeg =  Ap[node + 1] - Ap[node] ;
        }
    }

    /* Compute alpha, gamma, and the sum of all costs */
    if (options->c_is_all_ones == TRUE )
    {
        alpha = 1 ;
        gamma = 1 ;
        total_graph_weight = num_nodes ;
    }
    else
    {
        c = problem->graph->c ;
        alpha = c[0] ;
        gamma = c[0] ;
        total_graph_weight = 0 ;

        for (node = 0 ; node < num_nodes ; node++) {
            c_node = c [node] ;
            total_graph_weight += c_node ;
            if ( c_node > gamma ) gamma = c_node ;
            if ( c_node < alpha ) alpha = c_node ;
        }
    }
    problem->graph->min_vtx_cost = alpha ;
    problem->graph->max_vtx_cost = gamma ;
    problem->graph->total_vtx_cost = total_graph_weight ;


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




solve_result *VSP_Solve
(
    VSPINT             n,
    VSPLONG          *Ap,
    VSPINT           *Ai,
    VSPINT           *Ae, // Temporary argument while I'm refactoring
    Options     *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    Problem *probFinest, *probCoarsest, *P ;

    probFinest = Construct_Problem(n, Ap[n], Ai, Ap, Ae, options) ;
    set_default_problem_data(probFinest) ;

    /* Preprocessing */
    check_for_complete_graph(probFinest) ;
    preprocess_finest_problem(probFinest, options) ;

    #ifdef DEBUG
        check_problem(probFinest, options) ;
    #endif

    if(probFinest->graph->nz == 0)
    {
        solve_result *ret = Solve_Null_Graph(probFinest) ;
        Destroy_Coarse_Problems(probFinest) ;
        return ret ;
    }

    /*--------------------------------------------------------------------*/
    // Coarsen
    /*--------------------------------------------------------------------*/
    probCoarsest = Problem_Recursive_Coarsen(probFinest, options) ;
    State *X_l = NULL, *Y_l = NULL, *X_lp1 = NULL , *Y_lp1 = NULL,
    *X_0 = NULL, *Y_0 = NULL ;

    /*------------------------------------------------------------------------*/
    // Solve the coarsest level
    /*------------------------------------------------------------------------*/
    P = probCoarsest ;
    X_lp1 = Construct_State(P->graph->n) ;
    Y_lp1 = Construct_State(P->graph->n) ;
    set_initial_guess(X_lp1, Y_lp1, P, options) ;
    Solve_Level(X_lp1, Y_lp1, P, options) ;

    // Invariant: X_l, Y_l holds the solution at level l
    for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
    {
        X_l = Construct_State(P->graph->n) ;
        Y_l = Construct_State(P->graph->n) ;
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
    relaxed_to_vertex_separator(X_0, Y_0, probFinest) ;
    solve_result *result = Package_Result(X_0, Y_0, probFinest) ;
    Destroy_State(X_0) ;
    Destroy_State(Y_0) ;
    Destroy_Coarse_Problems(probFinest) ;
    return result ;
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
    bool dec = restart_decide(options->restart_decision, 0, 0, problem);
    if (!dec || problem->child_problem == NULL) return;

    Try_And_Keep_If_Better(X, Y, problem, strategy_centroid_restart, options) ;
}
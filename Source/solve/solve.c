#include "solve_internal.h"
#include "VSPQP.h"
#include "core.h"
#include "coarsen.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <threads.h>
#include "profiler.h"
#include "hubs.h"
#include "data_structures.h"

const char *G_Current_Problem_Name ;

bool restart_decide
(
    restart_decision           *dec,
    VSPFLOAT                  objIn, /* If not using rule, this is ignored */
    VSPFLOAT                 objOut, /* If not using rule, this is ignored */
    Problem                *problem  /* If using levels, this is ignored   */
) ;

/* Determine if we should resolve, and if so, do it */
static void Consider_Restart
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options
)
{
    VSPFLOAT ratio ;
    bool success ;
    bool dec = restart_decide(options->restart_decision, 0, 0, problem);
    if (!dec || problem->child_problem == NULL) return;

    void (*strategy)(State*, State*, Problem*, Options*) = NULL;
    switch (options->initial_guess)
    {
        case centroid_initial_guess:
            strategy = strategy_centroid_restart;
            break;
        case bisection_initial_guess:
            strategy = strategy_bisection_restart;
            break;
        case random_initial_guess:
            strategy = strategy_random_restart;
            break;
        default:
            fprintf(stderr, "ERROR: Unsupported starting guess for restart\n");
            return;
    }

    success = Try_And_Keep_If_Better(X, Y, problem, strategy, options, &ratio) ;

    // Record if restart was a success
    problem->stats.restart_success   = success ? 1 : 0 ;
    problem->stats.restart_obj_ratio = ratio ;
}

void Print_Stats
(
    Problem *Pfinest    
)
{
    FILE    *fp ;
    Problem  *P = Pfinest ;

    static bool first_time = true ;
    if(first_time)
    {
        first_time = false ;
        fp = fopen("output/stats.csv", "w") ;
        fprintf(fp,"id,lvl,n,edges,restart,ratio-sep\n") ;
    }
    else
    {
        fp = fopen("output/stats.csv", "a") ;
    }
    if(!fp) printf("Error: Could not open output/stats.csv") ;

    char *success = "success" ;
    char *fail    = "fail" ;
    char *na      = "" ;
    char *status ;

    while ( P != NULL)
    {
        if(P->stats.restart_success == EMPTY) status = na ;
        else status = (P->stats.restart_success) ? success : fail ;
        fprintf
        (
            fp, "%s,%lld,%lld,%lld,",
            G_Current_Problem_Name,
            PRI(P->level),
            PRI(P->graph->n),
            PRI(P->graph->nz / 2)
        ) ;

        if(status != na)
        {
            fprintf(fp, "%s,%f\n", status, P->stats.restart_obj_ratio) ;
        }
        else
        {
            fprintf(fp, ",\n") ;
        }

        P = P->child_problem ;
    }
    fclose(fp) ;
}

void Destroy_All_Problems(Problem *problemFinest)
{
    Problem *cp, *p = problemFinest ;
    while(p != NULL)
    {
        cp = p->child_problem ;
        Destroy_Problem(p) ;
        p = cp ;
    }
}


// Solve the VSP in the case that the graph has no edges
// Input: The graph has no edges, all vertex costs are 1
// Output: 
solve_result *Solve_Null_Graph
(
    VSPINT n,
    VSPINT u
)
{
    VSPINT i ;
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
    return result ;
}


solve_result *Package_Result(State *X, State *Y)
{
    VSPINT n = X->point->n ;
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
static void set_default_problem_data
(
    Problem* problem,
    Options *opts,
    VSPINT n_input
)
{
    VSPLONG i ;
    VSPINT   n = problem->graph->n ;
    VSPLONG nz = problem->graph->nz ;

    /* Set cost of each vertex to 1 */
    VSPINT* c = problem->graph->c ;
    for(i = 0 ; i < n ; i++) c[i] = 1 ;

    /* Set the upper balance constraint to 60%, and lower constraint to 1 */
    problem->u        = n_input * opts->u_ratio_finest ;
    problem->ubFinest = problem->u ;
    return ;
}

/* Compute alpha, gamma, sum of all costs, and set the is_original flag to
   true */
void preprocess_finest_problem
(
    Problem* problem,
    Options* options
)
{
    VSPINT node, num_nodes, c_node, alpha, gamma, total_graph_weight ;
    VSPINT* c ;
    VSPFLOAT ratio_cost_gamma ; // Ratio of notes with cost gamma
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
    problem->graph->max_deg = maxDeg ;

    /* Compute alpha, gamma, and the sum of all costs */
    if (problem->graph->c == NULL)
    {
        alpha = 1 ;
        gamma = 1 ;
        total_graph_weight = num_nodes ;
        ratio_cost_gamma = 1 ;
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

        VSPINT num_nodes_gamma = 0 ;
        for(node = 0 ; node < num_nodes ; node++)
        {
            if(c[node] == gamma) num_nodes_gamma++ ;
        }
        ratio_cost_gamma = ( (VSPFLOAT) (num_nodes_gamma) ) / num_nodes ;
    }
    problem->graph->min_vtx_cost = alpha ;
    problem->graph->max_vtx_cost = gamma ;
    problem->graph->total_vtx_cost = total_graph_weight ;
    problem->graph->ratio_max_cost = ratio_cost_gamma ;
    return ;
}


void check_for_complete_graph(VSPINT n, VSPLONG nz)
{
    long long upper_bnd_nz = (long long) n * n - n ;
    if( nz == upper_bnd_nz )
    {
        fprintf
        (
            stderr, "ERROR: Graph %s is complete\n", G_Current_Problem_Name
        ) ;
        exit(1) ;
    }
    else if(nz > upper_bnd_nz)
    {
        fprintf(stderr, "ERROR: nz > n^2 - n\n") ;
        exit(1) ;
    }
    return ;
}

// Attempt to move the root of the source heap (source = 0, 1, or 2) to the
// to state (to = 0 or 1)
static bool Move_Singleton
(
    int to, int source, State *states[2], min_heap_v2 *heaps[3], VSPINT u
)
{
    VSPINT cval, i ;
 
    // Check for empty heap
    MinHeap_See_Top_Key_v2(heaps[source], &cval) ;
    if(cval == EMPTY) return false ;

    // Check for infeasibility
    if(states[to]->weight + cval > u) return false ;

    // If transferring between shores, check if new move improces the imbalance
    if(source < 2)
    {
        VSPINT wto  = states[to]->weight;
        VSPINT woth = states[1 - to]->weight;    // the other shore
        VSPINT old_imb = ABS_VAL(wto - woth) ;
        VSPINT new_imb = ABS_VAL((wto + cval) - (woth - cval)) ;
        if (new_imb >= old_imb) return false ; // don't take a non-improving move
    }

    // Insert i into to
    MinHeap_Get_Top_Val_v2(heaps[source], &i) ;
    states[to]->point->int_part[i] = 1 ;
    states[to]->Az->int_part[i]++ ;
    states[to]->weight += cval ;

    if(source < 2)
    {
        // Remove i from other
        states[source]->point->int_part[i] = 0 ;
        states[source]->Az->int_part[i]-- ;
        states[source]->weight -= cval ;
    }

    MinHeap_Insert_v2(heaps[to], i) ;

    return true ;
}

// Note: ignores fractional components, for an easier implementation
// An iteration of the algorithm consts of moving a singleton
// from the separator or heavier shore (separator preferred)
// The iterations stop when one of the following stopping conditions are met:
// Stopping conditions:
// 1) Difference in weight is at most 1
// 2) No singletons in separator or heavier variable
// 3) Each singleton has moved N number of times (N = something like 2)
// 4) Moving any singleton from separator or heavier shore to lighter one is
//    infeasible
//
// Implementation detail: the third criterion is implicit. After a singleton
// has moved its maximum number of times, it is simply not added to the heap
void static Balance_Using_Singletons(State *X, State *Y, Problem* P)
{
    SET_FUNCTION_COUNTER ;
#ifdef DEBUG

    check_states(X, Y, P) ;
#endif
    /* ====================================================================== */
    // Variable declarations
    /* ====================================================================== */
    VSPINT i, j, to, source ;
    VSPINT  u = P->u ;
    VSPINT *c = P->graph->c ;
    VSPINT  n = P->graph->n ;
    VSPLONG *Ap = P->graph->Ap ;
    State *states[2] = {X, Y} ;  
    bool stop, moved ; 

    /* ====================================================================== */
    // Quit early if partition is already near-perfectly balanced
    /* ====================================================================== */
    if( imbalance(X, Y) <= 1 ) return ;

    /* ====================================================================== */
    // Build the 3 min heaps. 0 is for x, 1 for y, 2 for separator
    /* ====================================================================== */
    min_heap_v2 *heaps[3] ;
    for(i = 0 ; i < 3 ; i++) heaps[i] = MinHeap_Init_v2(c, n) ;

    // Insert singletons into the x heap, y heap, and singleton heap
    for(i = 0 ; i < n ; i++)
    {
        // Skip non-singletons
        if(Ap[i+1] - Ap[i] > 0) continue ;

        // Check if singleton is in intersection. Really shouldn't happen
        if(X->point->int_part[i] + Y->point->int_part[i] >= 2)
        {
            printf("WARNING: Singleton in intersection\n") ;

            // Remove the singleton from the bigger shore
            to = (states[0]->weight > states[1]->weight) ? 0 : 1 ;
            states[to]->point->int_part[i] = 0 ;
            states[to]->weight -= c[i] ;
        }

        // Skip singletons that have a fractional value
        if(i == X->point->ind_frac || i == Y->point->ind_frac) continue ;

        // Place i into the appropriate heap
        to = 2 ;
        for(j = 0 ; j < 2 ; j++)
        {
            if( states[j]->point->int_part[i] ) to = j ;
        }
        MinHeap_Insert_v2(heaps[to], i) ;
    }

    /* ====================================================================== */
    // Quit early if no singleton in separator or heavier shore
    /* ====================================================================== */
    // to is the bigger shore
    to = (states[0]->weight < states[1]->weight ? 1 : 0) ;
    if(heaps[2]->numNodes == 0 && heaps[to]->numNodes == 0) goto CLEANUP ;

    /* ====================================================================== */
    // Perform iterations
    /* ====================================================================== */
    stop = false ;

int loop = 0 ;
    while( !stop )
    {
        loop++ ;
        // Identify the smallest state
        to    = (states[0]->weight < states[1]->weight) ? 0 : 1 ;

        // Move a singleton from separator if possible
        moved = Move_Singleton(to, 2, states, heaps, u) ;
        
        // If no move from separator, attempt to move from the other shore
        if(!moved)
        {
            source = 1 - to ;
            moved = Move_Singleton(to, source, states, heaps, u) ;
        }

        // Stop if there was no move, or if balance is nearly perfect
        if(!moved || imbalance(states[0], states[1]) <= 1) stop = true ;
#ifdef DEBUG
    check_states(X, Y, P) ;
#endif
    }
    
    /* ====================================================================== */
    // Free memory
    /* ====================================================================== */

CLEANUP:
    for(i = 0 ; i < 3 ; i++) MinHeap_Free_v2(heaps[i]) ;

#ifdef DEBUG
    check_states(X, Y, P) ;
#endif
}


// Input: an almost binary point. Initial point at this level
// Output: all our improvements algorithms have been applied
void Solve_Level
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    static int count = 0 ;
    count++ ;

#ifdef DEBUG
    check_problem(problem, options) ;
#endif

    TIMER_START(T_TOTAL) ;
    
    Optimize(X, Y, problem, options) ;
    VSPFLOAT imbal_ratio = (VSPFLOAT) abs(X->weight - Y->weight) / problem->graph->total_vtx_cost ;

    if(   options->attempt_second_start
       && imbal_ratio > options->min_imbal_ratio_second_start
       && problem->level >= options->min_level_second_start)
    {
        Try_And_Keep_If_Better(X, Y, problem, strategy_explore_bisection, options, NULL) ;
    }

    Balance_Using_Singletons(X, Y, problem) ;

    TIMER_STOP(T_TOTAL) ;
    Timer_Flush(problem->level) ;
}

// Assumed: coarse X and Y have no fractional entries
void Uncoarsen_Point_NoSep
(
    Problem    *fine_graph,
    State          *fine_X,
    State          *fine_Y,
    State        *coarse_X,
    State        *coarse_Y
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT  fn       = fine_graph->graph->n ;
    VSPINT *matchmap = fine_graph->matching->matchmap ;

    VSPINT *fxInt = fine_X  ->point->int_part ;
    VSPINT *fyInt = fine_Y  ->point->int_part ;
    VSPINT *cxInt = coarse_X->point->int_part ;
    VSPINT *cyInt = coarse_Y->point->int_part ;


    for(VSPINT fv = 0 ; fv < fn ; fv++)
    {
        VSPINT cv = matchmap[fv] ;
        fxInt[fv] = cxInt[cv] ;
        fyInt[fv] = cyInt[cv] ;
    }

    fine_X->weight = coarse_X->weight ;
    fine_Y->weight = coarse_Y->weight ;

    Destroy_State(coarse_X) ;
    Destroy_State(coarse_Y) ;
}

solve_result *VSP_Solve
(
    VSPINT    n_input,
    VSPLONG  *Ap_input,
    VSPINT   *Ai_input,
    Options  *options,
    char     *name       // Optional. Problem name
)
{
    SET_FUNCTION_COUNTER ;
    G_Current_Problem_Name = name ;


    /*--------------------------------------------------------------------*/
    // Check for trivial inputs (no edges, or complete graphs)
    /*--------------------------------------------------------------------*/
    // Short-circuit if input graph is complete
    check_for_complete_graph(n_input, Ap_input[n_input]) ;

    // Short-circuit if input problem has no edges
    if(Ap_input[n_input] == 0)
    {
        return Solve_Null_Graph (n_input, n_input * options->u_ratio_finest) ;
    }

    /*--------------------------------------------------------------------*/
    // Remove hubs
    /*--------------------------------------------------------------------*/
    VSPINT   num_hubs ;
    VSPINT  *hubs ;
    VSPINT   n ;
    VSPINT  *Ai ;
    VSPLONG *Ap ;
    
    // Get hubs
    Get_Hubs(n_input, Ap_input, options->hub_treshold, &hubs, &num_hubs) ;

    // Every hub goes into the separator, so they only help if they're a small
    // minority. If there's a significant fraction of them, skip removal and
    // solve the original graph.
    bool use_hubs = (num_hubs > 0) && (num_hubs <= 0.3 * n_input) ;

    if(use_hubs)
    {
        Remove_Hubs(n_input, Ai_input, Ap_input, hubs, num_hubs, &n, &Ai, &Ap) ;
    }
    else
    {
        Ai = Ai_input ;
        Ap = Ap_input ;
        n  = n_input ;
        num_hubs = 0 ;
    }

    if(Ap[n] == 0)
    {
        solve_result *result = malloc(sizeof(*result)) ;
        result->x = calloc(n_input, sizeof(VSPINT)) ;
        result->y = calloc(n_input, sizeof(VSPINT)) ;

        VSPINT u = n_input * options->u_ratio_finest ;
        VSPINT x_count, y_count ;
        if(u < n / 2) { x_count = u ;     y_count = u ;         }
        else          { x_count = n / 2 ; y_count = n - n / 2 ; }

        VSPINT h = 0, placed = 0 ;
        for(VSPINT v = 0 ; v < n_input ; v++)
        {
            if(h < num_hubs && hubs[h] == v) { h++ ; continue ; }
            if      (placed < x_count)           result->x[v] = 1 ;
            else if (placed < x_count + y_count) result->y[v] = 1 ;
            placed++ ;
        }
        result->weight_x = x_count ;
        result->weight_y = y_count ;
        result->obj_val  = x_count + y_count ;

        free(hubs) ;
        free(Ai) ;
        free(Ap) ;
        return result ;
    }

    /*--------------------------------------------------------------------*/
    // Load problem into Problem struct
    /*--------------------------------------------------------------------*/
    bool have_global_optimum = false ;
    Problem *probFinest, *probCoarsest, *P ;

    probFinest = Construct_Problem(n, Ap[n], Ai, Ap, false, options) ;
    set_default_problem_data(probFinest, options, n_input) ;

    /* Preprocessing */
    preprocess_finest_problem(probFinest, options) ;

    #ifdef DEBUG
        check_problem(probFinest, options) ;
    #endif

    /*--------------------------------------------------------------------*/
    // Coarsen
    /*--------------------------------------------------------------------*/
    VSPINT *sol_sings ;
    VSPINT  n_sings ;
    probCoarsest = Problem_Recursive_Coarsen
                   (
                       probFinest, options, &sol_sings, &n_sings
                   ) ;
    State *X_l = NULL, *Y_l = NULL, *X_lp1 = NULL , *Y_lp1 = NULL,
    *X_0 = NULL, *Y_0 = NULL ;

    /*------------------------------------------------------------------------*/
    // Solve the coarsest level (or seed it from the singleton short-circuit)
    /*------------------------------------------------------------------------*/
    P = probCoarsest ;
    X_lp1 = Construct_State(P->graph->n) ;
    Y_lp1 = Construct_State(P->graph->n) ;

    if(sol_sings != NULL)
    {
        have_global_optimum = true ;

        VSPINT *c        = P->graph->c ;
        VSPINT  n_coarse = P->graph->n ;
        VSPINT  weight_x = 0 ;

        for(VSPINT i = 0 ; i < n_coarse ; i++)
        {
            X_lp1->point->int_part[i] = 0 ;
            Y_lp1->point->int_part[i] = 1 ;
        }
        for(VSPINT k = 0 ; k < n_sings ; k++)
        {
            VSPINT v = sol_sings[k] ;
            X_lp1->point->int_part[v] = 1 ;
            Y_lp1->point->int_part[v] = 0 ;
            weight_x += c[v] ;
        }
        X_lp1->weight = weight_x ;
        Y_lp1->weight = P->graph->total_vtx_cost - weight_x ;
    }
    else
    {
        set_initial_guess(X_lp1, Y_lp1, P, options) ;
        have_global_optimum =
                (X_lp1->weight + Y_lp1->weight == P->graph->total_vtx_cost) ;

        // Solve coarsest level if initial solution is not globally optimal
        if(!have_global_optimum) Solve_Level(X_lp1, Y_lp1, P, options) ;
    }

    // Invariant: X_l, Y_l holds the solution at level l.
    // On the singleton path the coarse solution has empty separator, so we
    // use the cheaper Uncoarsen_Point_NoSep (no Az computation).
    for(P = P->parent_problem ; P != NULL ; P = P->parent_problem)
    {

        X_l = Construct_State(P->graph->n) ;
        Y_l = Construct_State(P->graph->n) ;

        if(have_global_optimum)
        {
            Uncoarsen_Point_NoSep(P, X_l, Y_l, X_lp1, Y_lp1) ;
        }
        else
        {
            Uncoarsen_Point(P, X_l, Y_l, X_lp1, Y_lp1) ;
            Solve_Level(X_l, Y_l, P, options) ;
        }

        X_lp1 = X_l ;
        Y_lp1 = Y_l ;
    }

    // X_0, Y_0 is the solution at the finest, uncompressed level
    X_0 = X_lp1 ;
    Y_0 = Y_lp1 ;

    if(num_hubs)
    {
        Restore_Hubs(X_0, hubs, num_hubs) ;
        Restore_Hubs(Y_0, hubs, num_hubs) ;
        Place_Hubs_In_Shores
        (
            X_0, Y_0,
            hubs, num_hubs,
            Ai_input, Ap_input,
            probFinest->u
        ) ;
    }

    solve_result *result = Package_Result(X_0, Y_0) ;
    Destroy_State(X_0) ;
    Destroy_State(Y_0) ;
    free(sol_sings) ;

    Print_Stats(probFinest) ;
    Destroy_All_Problems(probFinest) ;

    free(hubs) ;
    if(num_hubs)
    {
        free(Ai) ;
        free(Ap) ;
    }

    return result ;
}
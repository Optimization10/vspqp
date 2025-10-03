#include "shared_api.h"

void VSP_set_default_options(VSP_options* options)
{
    /* Coarsening */
    options->use_coarsening             = true ;
    options->max_size_of_coarsest_graph = 10;
    options->do_community_match         = true ;
    options->matching_strategy          = HEMSR ;

    options->coarsenDec = NULL ;

    /* Initial guess */
    options->initial_guess                  = VSP_centroid_initial_guess;
    options->random_point_weight_proportion = 0.58 ;
    options->centroidVal = 0.6 ;

    /* Gamma reduction */
    options->use_gamma_reduction        = true ;
    options->num_gr_for_large_graphs    = 15 ;
    options->num_gr_for_small_graphs    = 30 ;
    options->max_num_vert_small_graph   = 150 ;

    options->sep_at_each_level          = false ;
    options->c_is_all_ones              = FALSE ;
    options->eps                        = 1E-8 ;


    /* Restart optimization at the current level */
    options->restart_decision = NULL ;

    /* Relaxing or constricting upper bound */
    options->relaxUB = false ;
    options->ratioUBLimit = 0.6 ;

    /* FM Options */
    options->FM_Limit = 400 ;

    /* ====================================================================== */
    // Mountain Climbing Options
    /* ====================================================================== */
    options->tol_LP_opt_checker = 1E-8 ;

    /* ====================================================================== */
    // Simple Coarsening options
    /* ====================================================================== */
    options->max_density     = 0.8 ;
    options->limit_num_nodes = 100 ;

    return ;
}

/* Set the cost of each vector to 1, each edge weight to 1, and
   the balance constraints as 60% of the total number of vertices             */
void VSP_set_default_problem_data(VSP_problem* problem) {
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT nz = problem->Graph->nz ;

    /* Set cost of each vertex to 1 */
    VSPINT* c = problem->Graph->c ;
    for(i = 0 ; i < n ; i++) c[i] = 1 ;
    problem->ratioMaxCostNodes    = 1 ;

    /* Set the edge weights to 1 */
    VSPINT* Ae = problem->Graph->Ae ;
    for(i = 0 ; i < nz ; i++) Ae[i] = 1 ;

    /* Set the upper balance constraint to 60%, and lower constraint to 1 */
    problem->u = floor(problem->Graph->n * 0.6);

    problem->ubFinest = floor(problem->Graph->n * 0.6) ;

    return ;
}

/* Compute alpha, gamma, sum of all costs, and set the is_original flag to
   true */
void    VSP_preprocess_finest_problem
(
    VSP_problem* problem,
    VSP_options* options
)
{
    VSPINT node, num_nodes, c_node, alpha, gamma, total_graph_weight ;
    VSPINT* c ;
    num_nodes = problem->Graph->n ;
    problem->was_coarsened_w_comm = false ;

    /* Finest problem has no parent. Lets us know when to stop refining */
    problem->parent_problem = NULL ;
    problem->level = 0 ;

    /* Compute the maximum degree */
    VSPINT maxDeg = 0 ;
    VSPINT *Ap    = problem->Graph->Ap ;
    for(node = 0 ; node < num_nodes ; node++)
    {
        if( Ap[node + 1] - Ap[node] > maxDeg )
        {
            maxDeg =  Ap[node + 1] - Ap[node] ;
        }
    }
    problem->maxDegree = maxDeg ;

    /* Compute alpha, gamma, and the sum of all costs */
    if (options->c_is_all_ones == TRUE )
    {
        alpha = 1 ;
        gamma = 1 ;
        total_graph_weight = num_nodes ;
    }
    else
    {
        c = problem->Graph->c ;
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
    problem->Graph->min_vtx_cost = alpha ;
    problem->Graph->max_vtx_cost = gamma ;
    problem->Graph->total_vtx_cost = total_graph_weight ;


    return ;
}
#include "core.h"
#include "coarsen_internal.h"
#include "coarsen.h"
#include <string.h>
#include <time.h>
#include <utils.h>
#include "profiler.h"

/* Coarsening works as follows: iterate through each vertex k of the coarse*/
/* Graph-> Let A and B be the vertices in the finer graph which*/
/* map to k. Iterate through all the neighbors of A. These will be neighbors*/
/* k, so add them to cAi. Munch is a variable which marches from 0 to the end */
/* of cAi. Increment munch as you go through the above-mentioned neighbors. */
/* Of course, at some point, you will come across B as one of the neighbors. */
/* Ignore it (do not increment munch). Now, go through the neighbors of */
/* B. It is possible that "neighbor" was a neighbor of vertex A, so we already*/ 
/* added it to cAi. We need to detect this. We do it via htable. htable[j] */
/* records the most recent index of cAi that j is at. So we can tell if j has */
/* already been recorded as a neighbor of k by looking if j is bigger than */
/* cAi[cAp[k]], which we recorded before as Ai.*/

/* Conditions to be met at the beginning of each loop: */
/* 1) munch = cAi[ cAp[k] ] (ps holds on to this value) */
/* 2) htable[j] = most recent index of cAi which contains j*/

/* Conditions to be met at the end of each loop: */
/* 1) munch = cAi[ cAp[k + 1] ]*/

static Graph* Graph_Coarsen
(
    Graph           *GraphFine,
    VSPINT            n_coarse,
    VSPINT           *matchmap,
    VSPINT        *invmatchmap,
    VSPINT           *matching,
    Options           *options
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;

        Check_Matrix
        (
            GraphFine->Ai,
            GraphFine->Ap,
            GraphFine->Ae,
            GraphFine->n,
            GraphFine->nz,
            TRUE,
            FALSE
        ) ;

    #endif

    VSPLONG *Ap_fine = GraphFine->Ap;
    VSPINT *Ai_fine  = GraphFine->Ai;
    VSPINT *Ae_fine  = GraphFine->Ae;
    VSPINT *c_fine   = GraphFine->c;


    /* Build the coarse graph */

    // Q: I am using the fine problem nz here, is this a problem?
    // A: No. The fine nz should be an upper bound on the coarse nz. And later
    //    on, we set problem_coarse->graph->nz = munch, which is the true coarse nz
    // Problem* problem_coarse = Construct_Problem(n_coarse, GraphFine->nz, options) ;
    Graph *Graph_Coarse = Construct_Graph(n_coarse, GraphFine->nz, NULL, NULL, NULL) ;

    VSPLONG *Ap_coarse = Graph_Coarse->Ap;
    VSPINT  *Ai_coarse = Graph_Coarse->Ai;
    VSPINT  *Ae_coarse = Graph_Coarse->Ae;
    VSPINT  *c_coarse  = Graph_Coarse->c;
    VSPINT   num_vtx_cost_gamma = 0 ;
    VSPLONG  munch = 0;

    /* Hashtable stores column pointer values. */
    VSPLONG* htable = malloc(n_coarse * sizeof(*htable));
    for (VSPINT i = 0; i < n_coarse; i++)
    {
        htable[i] = -1;
    }

    /* Coarse alpha and gamma, and maximum degree */
    VSPINT alpha_coarse = GraphFine->total_vtx_cost ;
    VSPINT gamma_coarse = 0;
    VSPINT maxDeg = 0 ;
    /* For each vertex in the coarse Graph-> */
    for (VSPINT k = 0; k < n_coarse; k++)
    {
        /* Load up the inverse matching */
        VSPINT v[3] = { EMPTY, EMPTY, EMPTY } ;
        v[0]     = invmatchmap[k] ;
        v[1]     = get_match(v[0], matching) ;
        if (v[0] == v[1])
        {
            v[1] = -1 ;
        }
        else
        {
            v[2] = get_match(v[1], matching);
            if (v[0] == v[2])
            {
                v[2] = -1;
            }
        }

        /* Right here, munch is the number of neighbors of 0, 1, 2, ..., k-1 */
        VSPLONG ps    = munch;
        Ap_coarse[k] = munch;

        /* Compute the degree of vertex k - 1 */
        if( k >= 1 )
        {
            if( Ap_coarse[k] - Ap_coarse[k-1] > maxDeg)
            {
                maxDeg = Ap_coarse[k] - Ap_coarse[k-1] ;
            }
        }

        VSPINT vertexCost = 0;
        for (VSPINT i = 0; i < 3 && v[i] != -1; i++)
        {
            /* Read the matched vertex and accumulate the vertex weight */
            /* and cost */
            VSPINT vertex = v[i];
            vertexCost += c_fine[vertex];
             

            for (VSPLONG p = Ap_fine[vertex]; p < Ap_fine[vertex + 1]; p++)
            {
                /* Inside of this double for loop, we are figuring out what
                   the neighbors of k in the coarse graph are. We do this by
                   looping through all neighbors of all preimages of k, and
                   adding them to the adjacency list of k.

                   We must ensure that we don't add the same coarse edge twice.
                   To do this, we keep htable. htable[r] is the most recent
                   position of r in the Ai array. If
                   htable[r] >= Ap_coarse[k], then the edge (k, r) has been
                   recorded already in the Ai segment for k. Otherwise, it has
                   not been recorded.

                   The value Ap[k] is saved as "ps"
                   "munch" just marches from 0 to nz_coarse. It is an indexor
                   for Ai_coarse.

                   */
                VSPINT neigh = Ai_fine[p] ;
                VSPINT toCoarsened = matchmap[neigh];
                if (toCoarsened == k)
                {
                    /* Delete self-edges */
                    continue;
                }

                /* Read the edge weight and accumulate the sum of edge weights*/
                /* technically it should be
                   Ae_fine[ Ai_fine[p] ];, but using the bugged
                   version Ae_fine[p] ]; gives both better time
                   and quality.
                   
                   UPDATE: No, the correct way is Ae_fine[p], not
                   Ae_fine[ Ai_fine[p] ]. The "argument" of
                   Ae is an index into the Ai array, which p is.
                   Ai[p] is the neighbor itself. But to find the edge weight
                   betweek i and j, you don't look up Ae[j]. Instead,
                   you go to the adjacency sublist in Ai corresponding to index
                   i (so the sublist from Ap[i] to Ap[i+1]), you locate j in
                   that sublist, you take note of the index of that j in the
                   Ap array (so the index being some p), and then you go to
                   index p in the Ae array
                    */
                VSPINT edgeWeight = Ae_fine[ p ] ;

                /* Check the hashtable before scattering. */
                VSPLONG cp = htable[toCoarsened];
                if (cp < ps) /* Hasn't been seen yet this column */
                {
                    htable[toCoarsened]        = munch;
                    Ai_coarse[munch]           = toCoarsened;
                    Ae_coarse[munch] = edgeWeight;
                    munch++;
                }
                /* If the entry already exists & we have edge weights,
                   sum the edge weights here. */
                else
                {
                    Ae_coarse[cp] += edgeWeight ;
                }
            }
        }
        c_coarse[k] = vertexCost;

        /* Update alpha, gamma, and num vertices with cost gamma */
        if (vertexCost > gamma_coarse)
        {
            gamma_coarse = vertexCost ;
            num_vtx_cost_gamma = 1 ;
        }
        else if (vertexCost == gamma_coarse) num_vtx_cost_gamma++ ;
        
        if (vertexCost < alpha_coarse)  alpha_coarse = vertexCost ;
    }
    Ap_coarse[n_coarse] = munch;

    /* Update max degree with the last vertex */
    if( Ap_coarse[n_coarse] - Ap_coarse[n_coarse - 1] > maxDeg)
    {
        maxDeg =  Ap_coarse[n_coarse] - Ap_coarse[n_coarse - 1] ;
    }

    /* Update coarse graph */
    Graph_Coarse->n              = n_coarse ;
    Graph_Coarse->nz             = munch ;
    Graph_Coarse->max_deg        = maxDeg ;
    Graph_Coarse->min_vtx_cost   = alpha_coarse ;
    Graph_Coarse->max_vtx_cost   = gamma_coarse ;
    Graph_Coarse->total_vtx_cost = GraphFine->total_vtx_cost ;
    Graph_Coarse->ratio_max_cost = ((VSPFLOAT) num_vtx_cost_gamma) / n_coarse ;

    /* Free memory */
    free(htable) ;

#ifdef DEBUG
    Check_Matrix
    (
        Graph_Coarse->Ai,
        Graph_Coarse->Ap,
        Graph_Coarse->Ae,
        Graph_Coarse->n,
        Graph_Coarse->nz,
        TRUE,   /* require_sorted -- now enforceable */
        FALSE
    ) ;
#endif

    #ifdef DEBUG_COARSENING

        VSPLONG end ;
        VSPINT i, v1, v2, k, r, s, is_edge ;
        VSPINT edge_weight_coarse_computed, edge_weight_coarse_claimed ;
        VSPLONG nz_coarse = Graph_Coarse->nz ;
        VSPLONG nz_fine   = GraphFine->nz ;
        VSPINT n_fine    = GraphFine->n ;


        /* Check that coarse gamma is correct */
        check_max_of_array(c_coarse, n_coarse, gamma_coarse) ;

        /* Check that Ap is (weakly) increasing */
        VSPINT Ap_increasing = TRUE ;
        for (i = 0 ; i < n_coarse ; i++) {
            if ( Ap_coarse [i] > Ap_coarse [i+1] ) Ap_increasing = FALSE ;
        }
        if (Ap_increasing == FALSE) {
            Error("Ap not increasing", DBG_INFO ) ;
        }

        /* Check that entries of Ai and Ap are within bounds */
        check_entries_are_between_bounds(Ai_coarse, nz_coarse, 0, n_coarse);

        /* If it's not too large, make the coarse graph explicit so it's
           easier for me to debug */
        VSPINT **A_fine = NULL, **A_coarse = NULL ;
        VSPINT **A_fine_weights = NULL, **A_coarse_weights = NULL ;
        if (n_fine < 5000) {
            A_fine = construct_explicit_matrix (
                                    NULL, Ai_fine, nz_fine, Ap_fine, n_fine) ;
            A_coarse = construct_explicit_matrix (
                            NULL, Ai_coarse, nz_coarse, Ap_coarse, n_coarse ) ;

            A_fine_weights = construct_explicit_matrix (
                        Ae_fine, Ai_fine, nz_fine, Ap_fine, n_fine) ;
            A_coarse_weights = construct_explicit_matrix (
              Ae_coarse, Ai_coarse, nz_coarse, Ap_coarse, n_coarse ) ;
        }

        /* === Check ==========================================================
           1) Check that whenever there is an edge (k, i) in the coarse graph,
           there is an edge between something in the preimage of k and
           something in the preimage of i.
           
           2) Check that the weight of (k, i) is the sum of the weights
           of all edges between the preimage of i and preimage of j

           3) Check that the cost and weight of every vertex in the coarse
              graph is the sum of the cost and weight (resp.) of the vertices
              in its preimage
           ================================================================== */

        for (k = 0 ; k < n_coarse ; k++) {
            end = Ap_coarse [k + 1] ;

            for (VSPLONG p = Ap_coarse [k] ; p < end ; p++ ) {

                /* There is an edge between k and i. Get the preimages
                   of k and i, i.e. the vertices in fine graph corresponding to
                   k and i */

                i = Ai_coarse [ p ] ;
                VSPINT* k_preimage = get_preimage_of_vertex(k, matching, invmatchmap) ;
                VSPINT* i_preimage = get_preimage_of_vertex(i, matching, invmatchmap) ;
                

                is_edge = FALSE ;

                /* Iterate through all pairs (v1, v2) where v1 is in the
                   preimage of k and v2 is in the preimage of i. Record whether
                   there is actually an edge between any of these pairs (there
                   should be, if there are no bugs in coarsening).
                   
                   Also, sum the edge weights of all existing edges of the form
                   (v1, v2). This should match the weight of edge (k, i)      */

                
                /* r, s are indices in the preimage arrays */
                v1 = k_preimage [ r = 0 ] ;
                edge_weight_coarse_computed = 0 ;
                VSPINT cost_k = 0 ;
                while ( v1 != EMPTY )
                {
                    cost_k += c_fine [v1] ;
                    v2 = i_preimage [ s = 0 ] ;
                    while (v2 != EMPTY)
                    {
                        if(there_is_an_edge(v1, v2, Ai_fine, Ap_fine))
                        {
                            is_edge = TRUE ;
                            edge_weight_coarse_computed += get_edge_weight(v1, v2, Ai_fine, Ap_fine, Ae_fine) ;
                        }
                        v2 = i_preimage [ ++s ] ;
                    }
                    v1 = k_preimage [ ++r ] ;
                }
                
                if (is_edge == FALSE)
                {
                    Error("Missing edge in coarse graph", DBG_INFO) ;
                }

                edge_weight_coarse_claimed = get_edge_weight(
                                                    k, i, Ai_coarse, Ap_coarse,
                                                    Ae_coarse) ;
                

                
                if ( edge_weight_coarse_computed != edge_weight_coarse_claimed)
                {
                    Error("Edge weight not correct.", DBG_INFO) ;
                }

                if (cost_k != c_coarse [k])
                {
                    Error("Cost of coarse vertex not correct", DBG_INFO) ;
                }

                /* Free preimages */
                free(k_preimage) ;
                free(i_preimage) ;
            }
        }

        /* === Check ========================================================
           1) The reverse direction of a check above. I.e. if there is an edge
              between any two members of any two preimages, then there is an
              edge connecting the corresponding coarse vertices.
           ================================================================== */
        VSPINT is_fine_edge, is_coarse_edge ;

        for (k = 0 ; k < n_fine ; k++) {
            VSPINT* matches_k = get_all_matches(k, matching) ;

            for (i = 0 ; i < n_fine ; i++) {

                VSPINT k_image = matchmap [k] ;
                VSPINT i_image = matchmap [i] ;

                if ( k_image == i_image ) continue ;
                else {
                    /* k, i map to different coarse vertices. Check if there
                       is an edge between the k LL and the i LL */
                    VSPINT* matches_i = get_all_matches(i, matching) ;

                    is_fine_edge =
                        check_for_edge_between_arrays_of_vertices (
                                matches_k, matches_i, 4, Ai_fine, Ap_fine ) ;

                    is_coarse_edge = there_is_an_edge (
                                    k_image, i_image, Ai_coarse, Ap_coarse) ;

                    if (is_fine_edge != is_coarse_edge) {
                        Error("Wrong edge in coarsening", DBG_INFO) ;
                    }

                    free (matches_i) ;
                }
            }
            free(matches_k) ;
        }


        if (A_fine != NULL) destruct_explicit_matrix(A_fine, n_fine) ;
        if (A_coarse != NULL) destruct_explicit_matrix(A_coarse, n_coarse) ;
        if (A_fine_weights != NULL) {
            destruct_explicit_matrix(A_fine_weights, n_fine) ;
        }
        if (A_coarse_weights != NULL) {
            destruct_explicit_matrix(A_coarse_weights, n_coarse) ;
        }
    #endif

    return Graph_Coarse;
}


// Generate a coarser version of the current problem
static Problem *Problem_Coarsen
(
    Problem *P,
    Options *options
)
{
    Problem *Q ;
    Problem *P_coarse ;
    Graph   *Graph_Coarse ;
    VSPINT   L ;

    Graph_Coarse = Graph_Coarsen
    (
        P->graph,
        P->matching->cn,
        P->matching->matchmap,
        P->matching->invmatchmap,
        P->matching->match, options
    ) ;

    P_coarse = Construct_Problem
    (
        Graph_Coarse->n,
        Graph_Coarse->nz,
        NULL, NULL, NULL,
        options
    ) ;

    // Free the graph coming from constructor
    Destroy_Graph(P_coarse->graph) ;
    P_coarse->graph = Graph_Coarse ;

    /* ====================================================================== */
    // Set other problem data
    /* ====================================================================== */
    P->child_problem         = P_coarse ;
    P_coarse->u              = P->u ;
    P_coarse->ubFinest       = P->ubFinest ;
    P_coarse->level          = P->level + 1 ;
    P_coarse->parent_problem = P ;

    // Set coarsest level num for each problem in the chain
    L = P_coarse->level ;
    Q = P_coarse ;
    while(Q != NULL)
    {
        Q->coarsest_level_num = L ;
        Q = Q->parent_problem ;
    }

    return P_coarse ;
}


static void SetUpperBound
(
    Problem *problem,
    Options *options
)
{
    if(options->relaxUB == false || problem->coarsest_level_num == 0 )
    {
        problem->u = problem->ubFinest ;
        return ;
    }

    VSPINT currLvl = problem->level ;
    VSPINT maxLvl  = problem->coarsest_level_num ;

    VSPFLOAT ratioUBLimit = options->ratioUBLimit ;

    // Linear Version
    VSPFLOAT gap          = ratioUBLimit - 0.6 ;
    VSPFLOAT ratioUb_Curr = 0.6 + currLvl * (gap / maxLvl) ;
    problem->u = (VSPINT) floor( ratioUb_Curr  * problem->graph->total_vtx_cost ) ;

    // Quadratic version
    // VSPFLOAT a  = (0.6 - ratioUBLimit) / (maxLvl * maxLvl) ;
    // VSPFLOAT ratioUb_Curr = a * (currLvl - maxLvl) * (currLvl - maxLvl) + ratioUBLimit ;
    // problem->u = (VSPINT) floor( ratioUb_Curr  * problem->graph->total_vtx_cost ) ;
    // printf("\nub: %f", problem->u / ((VSPFLOAT) problem->graph->total_vtx_cost )) ;

    return ;
}


// ORIGINAL VERSION
Problem *Problem_Recursive_Coarsen
(
    Problem *probFinest,
    Options    *options
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT    n_finest  = probFinest->graph->n ;
    VSPINT    l ;
    Problem  *P         = probFinest ;
    Problem  *coarse_P ;
    Problem  *coarsest_P ;

    l = 0 ; // Invariant: l = # of coarsenings
    while( 1 )
    {
        TIMER_START(T_COARSEN) ;
        // Match P
        Graph_Simple_Match
        (
            P->graph,
            P->matching,
            P->level,
            options
        ) ;

        // Coarsen P
        coarse_P = Problem_Coarsen(P, options) ;
 
        if(Graph_Density(coarse_P->graph->n, coarse_P->graph->nz) > 0.4
           || coarse_P->graph->n < 10)
        {
            coarsest_P = coarse_P ;
            break ;
        }
        else if(Is_Complete_Graph(coarse_P->graph->n, coarse_P->graph->nz))
        {
            Destroy_Problem(coarse_P) ;
            coarsest_P = P ;
            break ;
        }
        // Timer_Flush(G_Current_Problem_Name, l - 1) ;

        P = coarse_P ;
        l++ ;
        TIMER_STOP(T_COARSEN) ;
        Timer_Flush(l-1) ;
    }
    coarsest_P->child_problem = NULL ;

    /* Set upper bounds */
    P = coarsest_P ;
    while(P != NULL)
    {
       SetUpperBound(P, options) ;
       P = P->parent_problem ;
    }

    return coarsest_P ;
}
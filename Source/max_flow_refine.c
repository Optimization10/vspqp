#include "core.h"
#include <math.h>
#include "igraph.h"
#include "igraph_constants.h"
#include "igraph_flow.h"

#define OTHER_SIDE(i) (1 - (i))

#define INF_CAPACITY (gamma)
// We change the capacity of the infinite-capacity edges to gamma, because
// at the finest level, the network has unit capacities for every edge, and
// some algorithms run faster in this setting
// We can do this because this change to the capacities does not alter the
// set of feasible flows, (and hence does not alter the objective either).
// Consider the network with infinite edge capacities. Edges fall two four
// categories:
// 1) Edges whose endpoints are v_in for some vertex v of the original graph -
//    since any v_in has just one outgoing edge (to v_out) with capacity
//    <= gamma, no edge going into a v_in can carry more than gamma flow
// 2) Edges whose endpoints are t - these edges start at some v_out,
//    an the only incoming edge to v_out is some v_in. And edges of the form
//    (v_in, v_out) have capacity <= gamma. Hence, at most gamma flow
//    can come out of a v_out, so at most gamma flow can go along each edge
//    ending at t


// Gives 0, 1, or 2 if i is in side 0 (X side), side 1 (Y side), or separator
#define SIDE(i, sides) (2 - 2*(sides[0])[i] -  (sides[1])[i])


#define SS_INSERT(i, Adense, Asparse, numA)          \
    do { (Asparse)[i] = (numA) ;                     \
         (Adense)[(numA)++] = (i) ; } while (0)



                                            

static bool In_Sep(VSPINT i, VSPINT *sides[2])
{
    return sides[0][i] == 0 && sides[1][i] == 0 ;
}

static VSPLONG Degree(VSPINT i, VSPLONG *Ap)
{
    return Ap[i+1] - Ap[i] ;
}

// Return the number of vertices in A. A return value of 0 means that
// either the separator is empty, or its degree exceeds deg_budget
static VSPINT Get_A
(
    VSPINT  *Ai,
    VSPLONG *Ap,
    VSPINT   n,
    VSPINT   u,
    VSPINT  *c,
    VSPINT  *sides[2],
    VSPINT   wgts[3],     // assumed to be accurate
    VSPLONG  deg_budget,
    VSPINT  *Adense,      // can be garbage
    VSPINT  *Asparse      // assumed to be all -1's
)
{
    // Initialize variables
    VSPINT qstart = 0 ;
    VSPINT numA = 0 ;
    VSPINT side_budgets[2] = {u-wgts[2]-wgts[1],
                              u-wgts[2]-wgts[0]  } ;
    int side_enabled[3] = {TRUE, TRUE, TRUE} ;
    int deg_budget_remains = TRUE ;

    // Uninitialized
    int bi, bj ;
    VSPINT i, j, cost, deg ;
    VSPLONG p, end ;

    // Add S to A. Return early if separator exhausts all degree budget
    for(i = 0 ; i < n ; i++)
    {
        if(In_Sep(i, sides))
        {
            deg_budget -= Degree(i, Ap) ;
            if(deg_budget < 0)
            {
                // Clear Asparse and return
                for(j = 0 ; j < numA ; j++) Asparse[ Adense[j] ] = EMPTY ;
                return 0 ;
            }

            // Insert i into A and the queue
            Asparse[i] = numA ;
            Adense[numA++] = i ;
        }
    }

    // Do the BFS
    while( qstart < numA
           && (side_enabled[0] || side_enabled[1])
           && deg_budget_remains )
    {
        i = Adense[qstart++] ;
        bi = SIDE(i, sides) ;
        end = Ap[i+1] ;
        
        // If i is in the 0 or 1 side and that side is disabled, none of its
        // neighbors (which are in its own side or in the separator) can be
        // added to A. So skip i. i can however be in the separator, so
        // The side_enabled[2] is always TRUE, and this will make sure
        // we don't skip it, since sep can have neighbors on both sides
        if(side_enabled[bi] == FALSE) continue ;

        // Process neighbors of i
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;


            // Barrier 1: skip nodes of A (all separator nodes are in A already)
            if(Asparse[j] >= 0) continue ;

            bj = SIDE(j, sides) ; // i may be in the sep and j on a side

            // Permanently disable nodes entering from j's side if
            // the cost of j surpasses the remaining budget for that side
            cost = c[j] ;
            if(side_budgets[bj] - cost < 0) side_enabled[bj] = FALSE ;

            // Barrier 2: skip j, and all remaining neighbors of i, if j's
            // side became disabled
            if(side_enabled[bj] == FALSE)
            {
                if(bi < 2) break ;
                else      continue ;
            }

            // Barrier 3: stop building A if about to surpass the degree budget
            deg = Degree(j, Ap) ;
            if(deg_budget - deg < 0)
            {
                deg_budget_remains = FALSE ;
                break ;
            }

            // Insertion into A: finally, j passed all the barriers so it's
            // time to insert it in the queue, and decrease the budgets
            Asparse[j] = numA ;
            Adense[numA++] = j ;
            side_budgets[bj] -= cost ;
            deg_budget -= deg;
        }
    }
    return numA ;
}

/*
We are given a set of nodes A = (a_0, a_1, ..., a_{M-1}) 

In the flow network F, the nodes are labeled as follows:
s   = node 0
a_0 = nodes 1 (an in node) and 2 (an out node)
a_i = node 2*i+1 (an in node) and 2*i+2 (an out node)
a_{M-1} = node 2*M-1 (an in node) and 2*M (an out node)
t   = node 2*M + 1

Fix a node a_i in A, and let{a_j1, a_j2, ... a_jk} be its
neighbors in A. The arcs in F are:

2*i+1 (in version of a_i) --> 2*i+2 (out version of a_i), capacity c(a_i)
2*i+2 (out version of a_i) --> 2*j1 + 1 (in version of a_j1), infinite
s --> 2*i+1, infinite, if a_i has a neighbor in X outside of A
2*i+2 --> t, infinite, if a_i has a neighbor in Y outside of A

X \ A is contracted into s and Y \ A into t. This holds for every a_i,
whatever its side. The separator is in A, so a neighbor outside of A is
in X or in Y

*/

static void Generate_Max_Flow_Graph
(
    // Original graph and vertex costs
    VSPINT  *Ai,
    VSPLONG *Ap,
    VSPINT   n,
    VSPINT  *c,
    VSPINT   gamma,

    // The set A and the partition (X, S, Y)
    VSPINT  *Adense,
    VSPINT  *Asparse,
    VSPINT   numA,
    VSPINT  *sides[2],

    // Space to store the flow network in CSR
    // deg(A) = sum of the degrees of the vertices of A
    VSPINT  *Fi, // size at least 3 * numA + deg(A). Can be garbage
    VSPINT  *Fx, // size at least 3 * numA + deg(A). Can be garbage
    VSPLONG *Fp  // size at least 2 * numA + 3.      Can be garbage
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, j, k, u ;
    VSPINT M = 2 * numA + 2 ; // number of vertices in F 
    VSPLONG p, end ;
    VSPLONG t = 0 ; // walks along Fi and Fx
    int to_s, to_t ;

    // Initialize Fp. Fp[x+1] counts the arcs leaving node x of F. In-nodes
    // have just one arc leaving from them. For s and the out nodes, the
    // count is initialized to 0. No arc leaves t
    Fp[0] = 0 ;
    Fp[1] = 0 ;
    for(i = 1 ; i <= numA ; i++)
    {
        Fp[2*i]   = 1 ;  // in node 2i-1
        Fp[2*i+1] = 0 ;  // out node 2i
    }
    Fp[M] = 0 ;

    // Walk A to find the endpoints of arcs leaving s, and the initial points
    // of arcs entering t. Also count the arcs leaving each out node.
    // Only the sides of the neighbors outside of A matter, not the side of i
    for(k = 0 ; k < numA ; k++)
    {
        // i is a vertex in A. It corresponds to vertices 2k+1 and 2k+2 in F
        i    = Adense[k] ;
        to_s = FALSE ;
        to_t = FALSE ;

        end = Ap[i+1] ;
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;

            if(Asparse[j] != EMPTY)      Fp[2*k+3]++ ;  // arc (2k+2, in of j)
            else if(SIDE(j, sides) == 0) to_s = TRUE ;
            else                         to_t = TRUE ;  // j is in Y
        }

        // Make the arc (0, 2k+1) with infinite capacity. The arcs of s come
        // first in Fi, so they can be written now
        if(to_s)
        {
            Fp[1]++ ;
            Fi[t]   = 2*k+1 ;
            Fx[t++] = INF_CAPACITY ;
        }

        // The arc (2k+2, t) is written in the second walk
        if(to_t) Fp[2*k+3]++ ;
    }

    // Finalize Fp by accumulating its entries. Fp[x] is now the position in
    // Fi of the first arc leaving node x
    for(i = 1 ; i < M + 1 ; i++) Fp[i] += Fp[i-1] ;

    // Walk A, which is the same as walking the vertices of F in order from
    // least to greatest
    for(k = 0 ; k < numA ; k++)
    {
        to_t = FALSE ;
        // i is a vertex in A. It corresponds to vertices 2k+1 and 2k+2 in F
        i = Adense[k] ;

        // We are processing vertices 2k+1 and 2k+2 of F. First, make the arc
        // (2k+1, 2k+2) whose capacity is the cost of the underlying vertex of A
        Fi[t]   = 2*k + 2 ;
        Fx[t++] = c[i] ;
        
        // Walk the neighbors of i to obtain the destinations of arcs
        // from 2k + 2
        end = Ap[i+1] ;
        for(p = Ap[i] ; p < end ; p++)
        {
            j = Ai[p] ;

            u = Asparse[j] ; // u is position of j in Adense
            if( u != EMPTY ) 
            {
                // u corresponds to 2*u+1 and 2*u+2
                // Make the arc (2k+2, 2u+1)
                Fi[t]   = 2*u + 1 ;
                Fx[t++] = INF_CAPACITY ;
            }
            else if(SIDE(j, sides) == 1) to_t = TRUE ;
        }

        // Add the arc (2k+2, t) if i has a neighbor outside of A in Y. This
        // includes separator vertices
        if(to_t)
        {
            Fi[t]   = M - 1 ;
            Fx[t++] = INF_CAPACITY ;
        }
    }

    // The arcs written must match the counts in Fp
    VSP_ASSERT(t == Fp[M]) ;

    return ;
}

static void Get_Partition
(
    // The Network
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,  // number of vertices in F
    int     *part,
    VSPINT   seed, // start the BFS here   
    VSPINT   seed_part_num
)
{
        
}

static void Get_Max_Flow
(
    // The Network
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,  // number of vertices in F
    int     *part
)
{
    VSPINT i, j ;
    VSPLONG p, end ;
    VSPINT eF = Fp[vF] ; // number of edges in F

    // ========================================================================
    // Store the network edges and edge capacities in an igraph vector
    // ========================================================================
    igraph_vector_int_t edges ;
    igraph_vector_t capacities ;

    // Each edge is stored by its start and end points. So we need space for
    // twice the number of edges
    igraph_vector_int_init(&edges, 2*eF) ;
    igraph_vector_init(&capacities, eF) ; // one capacity per edge

    for(i = 0 ; i < vF ; i++)
    {
        end = Fp[i+1] ;
        for(p = Fp[i] ; p < end ; p++)
        {
            j = Fi[p] ;
            VECTOR(edges)[2*p]     = i ;
            VECTOR(edges)[2*p + 1] = j ;
            VECTOR(capacities)[p]  = Fx[p] ;
        }
    }

    // ========================================================================
    // Create an igraph network
    // ========================================================================
    igraph_t F ;
    igraph_create(&F, &edges, vF, IGRAPH_DIRECTED) ;


    // ========================================================================
    // Get the optimal flow
    // ========================================================================
    igraph_vector_t flow ;
    igraph_vector_init(&flow, eF) ; // one flow per edge
    igraph_maxflow
    (&F, NULL, &flow, NULL, NULL, NULL, 0, vF - 1,  &capacities, NULL) ;

    // ========================================================================
    // Get the the two extreme separators; one with the smallest s side,
    // and the other with the smallest t side
    // ========================================================================

}

void Max_Flow_Refine
(
    State   *X,
    State   *Y,
    Problem *problem
)
{
    VSPINT numA, M, i ;
    VSPINT  *Ai       = problem->graph->Ai ;
    VSPLONG *Ap       = problem->graph->Ap ;
    VSPINT   n        = problem->graph->n ;
    VSPINT  *sides[2] = {X->point->int_part, Y->point->int_part} ;
    VSPINT   wgts[3]  = {X->weight, Y->weight,
                         problem->graph->total_vtx_cost -X->weight -Y->weight} ;
    VSPINT   u        = problem->u ;
    VSPINT  *c        = problem->graph->c ;
    VSPINT   gamma    = problem->graph->max_vtx_cost ;
    
    VSPINT   deg_budget = 256 * sqrt(n) ;
    VSPINT  *Adense     = malloc(n * sizeof(*Adense)) ;
    VSPINT  *Asparse    = malloc(n * sizeof(*Asparse)) ;
    for(i = 0 ; i < n ; i++) Asparse[i] = EMPTY ;

    VSPINT  *Fi = malloc( (2 * Ap[n] + 2) * sizeof(*Fi)) ;
    VSPINT  *Fx = malloc( (2 * Ap[n] + 2) * sizeof(*Fx)) ;
    VSPLONG *Fp = malloc( (2 * n + 3) * sizeof(*Fp)) ;


    for(int num_its = 0 ; num_its < 3 ; num_its++)
    {
        // Get the area A. Quit if no area
        numA = Get_A(Ai, Ap, n, u, c, sides, wgts, deg_budget, Adense, Asparse);
        if(numA == 0) break ;
        

        // Get the flow graph corresponding to the area
        Generate_Max_Flow_Graph
        (Ai, Ap, n, c, gamma, Adense, Asparse, numA, sides, Fi, Fx, Fp) ;
        M = 2*numA + 2 ; // number of vertices in s-t flow graph

        // Construct the residual network of the final flow

        // Read off the two extreme min cuts from the residual 
    }

    // Free memory
    free(Fi) ;
    free(Fx) ;
    free(Fp) ;
    free(Adense) ;
    free(Asparse) ;
}
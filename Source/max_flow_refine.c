#include "core.h"
#include <math.h>
#include "igraph.h"
#include "igraph_constants.h"
#include "igraph_flow.h"
#include "igraph_vector.h"

#define OTHER_SIDE(i) (1 - (i))

#define INF_CAPACITY (gamma + 1)
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

// Return the number of vertices in R. A return value of 0 means that
// either the separator is empty, or its degree exceeds deg_budget
static VSPINT Get_R
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


// Given the area A of the original graph, generate the max flow graph
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



static void Get_Partition_Smallest_s
(
    // The Network
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,  // number of vertices in F

    // BFS info
    igraph_vector_t flow,

    // Partition of A
    int *part
)
{
    VSPINT i, j ;
    VSPLONG p, end ;
    igraph_real_t residual_flow ;
    VSPINT vA = (vF - 2) / 2 ; // number of vertices in A
    VSPINT *q       = malloc(vF * sizeof(*q)) ;
    q[0] = 0 ;
    VSPINT start_q  = 0 ;
    VSPINT end_q    = 0 ;
    int    *visited = calloc(vF, sizeof(*visited)) ;
    visited[0]   = 1 ;

    // If doing the BFS from s, initially place all vertices of A onto the
    // 1 side, which corresponds to t. If doing the BFS from t, initially place
    // all vertices of A onto the 0 side

    // Mark each node in partition 1 initially
    for(i = 0 ; i < vA ; i++) part[i] = 1 ;


    while(start_q <= end_q)
    {
        i = q[start_q++] ;

        // Map i back to a vertex j of A
        j = (i-1) / 2 ;
        part[j] = 0 ;

        // Visit the out-neighbors of i in the residual graph R
        
        // If i is a v_in node, then its only potential out-neighbors is v_out,
        // which is the same neighbor it has in F. This neighbor is Fi[Fp[i]]

        // If i is a v_out node, then in the residual graph its potential
        // out-neighbors are the same neighbors it has in F, plus the s node
        // (which, even if it is an out-neighbor in the residual graph, has been
        // visited so would be skipped anyway), plus v_in. v_in is the more
        // difficult casel; v_in is an out-neighbor of v_out iff the flow
        // across (v_in, v_out) is non-zero. So we need to know if this
        // flow is non-zero in order to know all out-neighbors of v_out

        
        // If i = s, then 

        // We don't have to worry about i = t because it will never be
        // reached


        // Place neighbors of i into the queue
        end = Fp[i+1] ;
        for(p = Fp[i] ; p < end ; p++)
        {
            j = Fi[p] ;

            // Compute residual flow on edge (i, j)
            // residual_flow
        }
    }


    free(q) ;
    free(visited) ;
}





// Given a flow netowrk, get an optimal flow
static igraph_vector_t Get_Max_Flow
(
    // The Network
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF   // number of vertices in F
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
    return flow ;
}

// Read off the two partitions P, C, Q of R based on the flow.
// Rpart1 contains the partition with smallest P. It is obtained by setting P
// to be all vertices reacheable from p in the residual graph
// Rpart2 contains the partition with smallest Q. The Q consists of all vertices
// from which q is reacheable in the residual graph
static void Apply_Flow_Refinement
(
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,   // numver of vertices
    VSPINT   gamma,
    igraph_vector_t flow,
    VSPINT *Rdense,
    VSPINT  *c // cost vector of the original graph
)
{
    VSPINT p, q, vR, i, k, M, j, v_in, v_out, u_in, u, v, f ;
    VSPLONG *Rp, *ptrs ;
    VSPLONG t, end ;
    VSPINT *Ri, *Rx ;
    p = 0 ;
    q = vF-1 ;

    // ========================================================================
    // Build the residual graph (TODO: allocate this memory only once, before)
    // ========================================================================
    vR = vF ;

    /* --- Build Rp array --------------------------------------------------- */
    // First, store the number of neighbors of vertex i in Rp[i+1]
    // Account for neighbors of p. For any arc (p, i), i gets another
    // outgoing arc in the residual graph
    Rp = calloc(vF, sizeof(*Rp)) ;
    for(t = Fp[0] ; t < Fp[1] ; t++) Rp[ Fi[t]+1 ]++ ;


    for(i = 1 ; i < vR - 1 ; i++)
    {
        // Even nodes not in {p, q} get an additional residual arc
        Rp[i+1] += Fp[i+1] - Fp[i] + (i % 2 == 0) ;

        // Count arcs going towards q
        for(t = Fp[i+1] ; t < Fp[i] ; t++)
        {
            if(Fi[t] == vR - 1) Rp[vF]++ ;
        }
    }
    for(i = 0 ; i < vR + 1; i++) Rp[i] += Rp[i-1] ; // accumulate Rp

    /* --- Build Ri and Rx arrays ------------------------------------------- */
    // ptrs[i] is initialized as the starting Ri-index of the adjacency list
    // of vertex i
    Ri   = malloc(Rp[vF] * sizeof(*Ri)) ;
    Rx   = malloc(Rx[vF] *  sizeof(*Rx)) ;
    ptrs = malloc(vR * sizeof(*ptrs)) ;
    for(i = 0 ; i < vR ; i++) ptrs[i] = Rp[i+1] ;

    /* --- (1) arcs of the form (p, v_in) ----------------------------------- */
    for(t = Fp[p] ; t < Fp[p+1] ; t++)
    {
        v_in = Fi[p] ; // t is the index in F of the edge (p, v_in)

        // Add the arc (p, v_in, infinity)
        Ri[ ptrs[p]   ] = v_in ;
        Rx[ ptrs[p]++ ] = INF_CAPACITY ;

        // Add the arc (v_in, p, flow(p, v_in))
        f = VECTOR(flow)[t] ;
        Ri[ ptrs[v_in]   ] = p ;
        Rx[ ptrs[v_in]++ ] = f ;
    }

    // Loop through the arcs of F
    for(i = 1 ; i < vF ; i++)
    {
        for(t = Fp[i] ; Fp[i+1] ; t++)
        {
            j = Fi[t] ;

            /* --- (2) if (i, j) = (v_in, v_out) ---------------------------- */
            if( (i % 2) && j == i + 1)
            {
                v_in  = i ;
                v_out = j ;
                
                // Add arc (v_in, v_out, c_v - flow(v_in, v_out))
                Ri[ptrs[v_in]  ] = v_out ;
                Rx[ptrs[v_in]++] = c[ Rdense[(i-1) / 2] ] - VECTOR(flow)[t]  ;

                // Add arc (v_out, v_in, flow(v_in, v_out)) 
                Ri[ptrs[v_out]  ] = v_in ;
                Rx[ptrs[v_out]++] = VECTOR(flow)[t] ;
                continue ;
            }

            /* --- (3) if (i, j) = (v_out, q) ------------------------------- */
            if(j == q)
            {
                v_out = j ;

                // Add arc (v_out, q, infty)
                Ri[ptrs[v_out]  ] = q ;
                Rx[ptrs[v_out]++] = INF_CAPACITY  ;

                // Add arc (q, v_out, f(v_out, q))
                Ri[ptrs[q]  ] = v_out ;
                Rx[ptrs[q]++] = VECTOR(flow)[t] ;
                continue ;
            }

            /* --- (4) if (i, j) = (v_out, u_in) ---------------------------- */
            v_out = i ;
            u_in = j ;

            // Add arc (v_out, u_in, infty)
            Ri[ptrs[v_out]  ] = u_in ;
            Rx[ptrs[v_out]++] = INF_CAPACITY  ;

            // Add arc (u_in, v_out, flow(v_out, u_in) )
            Ri[ptrs[u_in]  ] = v_out ;
            Rx[ptrs[u_in]++] =  VECTOR(flow)[t] ;
        }
    }

    // ========================================================================
    // Do a BFS from p to obtain cut1
    // ========================================================================
    int *cut1 = calloc(vF, sizeof(*cut1)) ;
    cut1[0] = 1 ;

    VSPINT res_flow, qstart, qend ;
    VSPINT *queue = malloc(vF * sizeof(queue)) ;
    
    // initialize the queue with p
    qstart        =  0 ;
    qend          = -1 ;
    queue[++qend] = p ;
    cut1[p]       = 1 ;
    while(qstart <= qend)
    {
        i = queue[qstart++] ;
        for(t = Rp[i] ; t < Rp[i+1] ; t++)
        {
            // Skip neighbors with no residual flow
            res_flow = Rx[t] ;
            if(!res_flow || cut1[j]) continue ;
            j = Ri[t] ;
            queue[++qend] = j ;
            cut1[j] = 1 ;
        }
    }

    // ========================================================================
    // Get the transpose of R in order to be able to follow backward arcs
    // ========================================================================
    // Every arc t = (i, j) of R becomes the arc (j, i) of RT, with the same
    // capacity Rx[t]. So i can reach j in R along t iff RT has (j, i) with
    // RTx > 0, and a BFS from q in RT finds all nodes that can reach q in R.
    //
    // The count of node j is stored two places ahead, in RTp[j+2]. After
    // accumulating, RTp[j+1] is the start of the list of j, and serves as its
    // write pointer. After the write pass, RTp[j+1] has moved to the end of
    // the list of j, which is the start of the list of j+1. So RTp[0..vF] is
    // the final pointer array of RT, with no separate ptrs array needed
    VSPLONG *RTp = calloc(vF + 2, sizeof(*RTp)) ;
    VSPINT  *RTi, *RTx ;

    // Make a pass to record the number of entries per row
    for(i = 0 ; i < vF ; i++)
    {
        for(t = Rp[i] ; t < Rp[i+1] ; t++)
        {
            j = Ri[t] ;
            RTp[j+2]++ ;
        }
    }

    // accumulate RTp. RTp[vF+1] is now the number of arcs of RT (= of R)
    for(i = 1 ; i < vF + 2 ; i++) RTp[i] += RTp[i-1] ;
    RTi = malloc(RTp[vF+1] * sizeof(*RTi)) ;
    RTx = malloc(RTp[vF+1] * sizeof(*RTx)) ;

    // Make another pass to write RTi and RTx
    for(i = 0 ; i < vF ; i++)
    {
        for(t = Rp[i] ; t < Rp[i+1] ; t++)
        {
            j = Ri[t] ;
            RTi[ RTp[j+1]   ] = i ;
            RTx[ RTp[j+1]++ ] = Rx[t] ;
        }
    }

    // ========================================================================
    // Do a BFS from q in RT to obtain cut2
    // ========================================================================
    // cut2[i] = 1 iff q is reachable from i in R, so 1 marks the q side here
    // (in cut1, 1 marks the p side). When the BFS stands at i and follows the
    // arc (i, j) of RT with RTx > 0, R has the arc (j, i) with positive
    // residual capacity, so j reaches i in R, and hence reaches q
    int *cut2 = calloc(vF, sizeof(*cut2)) ;

    // initialize the queue with q (the queue of the first BFS is reused)
    qstart        =  0 ;
    qend          = -1 ;
    queue[++qend] = q ;
    cut2[q]       = 1 ;
    while(qstart <= qend)
    {
        i = queue[qstart++] ;
        for(t = RTp[i] ; t < RTp[i+1] ; t++)
        {
            // Skip arcs with no residual flow, and nodes already visited
            j = RTi[t] ;
            res_flow = RTx[t] ;
            if(!res_flow || cut2[j]) continue ;
            queue[++qend] = j ;
            cut2[j] = 1 ;
        }
    }


    // ========================================================================
    // Figure out the cost of the shores given the separating partition P, C, Q
    // from cut 1
    // ========================================================================
    // Go through the cut. For a vertex v, if both v_in and v_out, 
    int *separation1 = malloc(vF * sizeof(*separation1)) ;
    int *separation2 = malloc(vF * sizeof(*separation2)) ;

    

    // Free memory
    free(separation1) ;
    free(separation2) ;
    free(cut2) ;
    free(RTp) ;
    free(RTi) ;
    free(RTx) ;
    free(queue) ;
    free(cut1) ;
    free(ptrs) ;
}

// Refine the partition represented by X and Y
void Max_Flow_Refine
(
    State   *X,
    State   *Y,
    Problem *problem
)
{
    igraph_vector_t flow ;
    VSPINT numA, vF, i ;
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

    int *Rpart1  = malloc(n * sizeof(*Rpart1)) ;
    int *Rpart2  = malloc(n * sizeof(*Rpart2)) ;


    for(int num_its = 0 ; num_its < 3 ; num_its++)
    {
        // Get the area A. Quit if no area
        numA = Get_R(Ai, Ap, n, u, c, sides, wgts, deg_budget, Adense, Asparse);
        if(numA == 0) break ;
        

        // Get the flow graph corresponding to the area
        Generate_Max_Flow_Graph
        (Ai, Ap, n, c, gamma, Adense, Asparse, numA, sides, Fi, Fx, Fp) ;
        vF = 2*numA + 2 ; // number of vertices in s-t flow graph

        // Get the max flow
        flow = Get_Max_Flow(Fi, Fx, Fp, vF) ;

        // Read off the partition of R off the residual network


        // Read off the two extreme min cuts from the residual 
    }

    // Free memory
    free(Fi) ;
    free(Fx) ;
    free(Fp) ;
    free(Adense) ;
    free(Asparse) ;
    free(Rpart1) ;
    free(Rpart2) ;
}
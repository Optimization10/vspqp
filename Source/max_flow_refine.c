#include "core.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include "igraph.h"
#include "igraph_constants.h"
#include "igraph_flow.h"
#include "igraph_vector.h"

#define OTHER_SIDE(i) (1 - (i))

// Rules that skip the max flow refinement when it cannot pay for itself. Both
// are checked at the start of every pass, on the cost C(S) of the separator
// (= |S| at the finest level, where all vertex costs are 1)
// 1) C(S) < MF_MIN_SEP: the max flow value is at most C(S), so at most a few
//    vertices can be saved
// 2) C(S) > MF_MAX_SEP_FACTOR * sqrt(n): the region A has a degree budget of
//    256 sqrt(n), and S alone takes up a large part of it, so A is barely
//    larger than S and the max flow can only make small local changes
#define MF_MIN_SEP        3
#define MF_MAX_SEP_FACTOR 10.0

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

    // Free the igraph objects. The caller frees flow
    igraph_destroy(&F) ;
    igraph_vector_destroy(&capacities) ;
    igraph_vector_int_destroy(&edges) ;
    return flow ;
}

// Read off the two partitions P, C, Q of R based on the flow.
// Rpart1 contains the partition with smallest P. It is obtained by setting P
// to be all vertices reacheable from p in the residual graph
// Rpart2 contains the partition with smallest Q. The Q consists of all vertices
// from which q is reacheable in the residual graph
// Returns 1 if a vertex moved (X, S, Y changed), 0 if not
static int Apply_Flow_Refinement
(
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,   // numver of vertices in F
    VSPINT   gamma,
    igraph_vector_t flow,
    VSPINT *Rdense,
    VSPINT  *c, // cost vector of the original graph
    VSPINT   cost_p, // cost of X \ R
    VSPINT   cost_q,  // cost of Y \ R
    State *  X,
    State *  Y,
    Problem *problem
)
{
    SET_FUNCTION_COUNTER ;

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
    // First, store the number of neighbors of vertex i in Rp[i+1]. Every arc
    // (i, j) of F gives two arcs of R: (i, j) leaving i, and (j, i) leaving j.
    // So the number of arcs leaving i in R is the number of arcs leaving i in
    // F plus the number of arcs entering i in F. This covers all four cases
    // below: (p, v_in), (v_in, v_out), (v_out, q), and (v_out, u_in)
    Rp = calloc(vR + 1, sizeof(*Rp)) ;
    for(i = 0 ; i < vR ; i++)
    {
        Rp[i+1] += Fp[i+1] - Fp[i] ;                       // arcs leaving i
        for(t = Fp[i] ; t < Fp[i+1] ; t++) Rp[ Fi[t]+1 ]++ ; // arcs entering Fi[t]
    }
    for(i = 1 ; i < vR + 1; i++) Rp[i] += Rp[i-1] ; // accumulate Rp

    /* --- Build Ri and Rx arrays ------------------------------------------- */
    // ptrs[i] is initialized as the starting Ri-index of the adjacency list
    // of vertex i
    Ri   = malloc(Rp[vR] * sizeof(*Ri)) ;
    Rx   = malloc(Rp[vR] * sizeof(*Rx)) ;
    ptrs = malloc(vR * sizeof(*ptrs)) ;
    for(i = 0 ; i < vR ; i++) ptrs[i] = Rp[i] ;

    /* --- (1) arcs of the form (p, v_in) ----------------------------------- */
    for(t = Fp[p] ; t < Fp[p+1] ; t++)
    {
        v_in = Fi[t] ; // t is the index in F of the edge (p, v_in)

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
        for(t = Fp[i] ; t < Fp[i+1] ; t++)
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
                v_out = i ;

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

    // Every list of R must be exactly full: the arcs written must match the
    // counts in Rp
    for(i = 0 ; i < vR ; i++) VSP_ASSERT(ptrs[i] == Rp[i+1]) ;

    // ========================================================================
    // Do a BFS from p to obtain cut1
    // ========================================================================
    // In cut1:
    // 1 -- P hat
    // 0 -- Q hat
    int *cut1 = calloc(vF, sizeof(*cut1)) ;
    cut1[0] = 1 ;

    VSPINT res_flow, qstart, qend ;
    VSPINT *queue = malloc(vF * sizeof(*queue)) ;
    
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
            // Skip arcs with no residual flow, and nodes already visited
            j = Ri[t] ;
            res_flow = Rx[t] ;
            if(!res_flow || cut1[j]) continue ;
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

    // In cut2:
    // 0 -- P hat
    // 1 -- Q hat
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
    // In PCQ_1 array, indexed by the position v of a vertex in R (Rdense):
    // 0 -- in P
    // 1 -- in Q
    // 2 -- in C
    VSPINT sizeR = (vF - 2) / 2 ;
    int *PCQ_1 = malloc(sizeR * sizeof(*PCQ_1)) ;
    VSPINT sep1Pcost = cost_p ;
    VSPINT sep1Qcost = cost_q ;
    VSPINT sep1Ccost = 0 ;

    for(i = 1 ; i <= sizeR ; i++)
    {
        v     = i-1 ;
        v_in  = 2*i - 1 ;
        v_out = v_in + 1 ;

        if( cut1[v_in] == 1 && cut1[v_out] == 1) // v goes in P
        {
            PCQ_1[v]   = 0 ;
            sep1Pcost += c[ Rdense[v] ] ;
        }
        else if(cut1[v_in] == 0 && cut1[v_out] == 0) // v goes in Q
        {
            PCQ_1[v]   = 1 ;
            sep1Qcost +=  c[ Rdense[v] ] ;
        }
        else // v goes in C
        {
            PCQ_1[v]   = 2 ;
            sep1Ccost +=  c[ Rdense[v] ] ;
        }
    }
    VSPINT sep1imb = ABS_VAL(sep1Pcost - sep1Qcost) ;


    // ========================================================================
    // Figure out the cost of the shores given the separating partition P, C, Q
    // from cut 2
    // ========================================================================
    // Here 1 marks the q side: v_in and v_out both reach q -> Q, neither
    // does -> P, only v_out does -> C (only v_in reaching q is impossible)
    //
    // In PCQ_2 array, indexed by the position v of a vertex in R (Rdense):
    // 0 -- in P
    // 1 -- in Q
    // 2 -- in C
    int *PCQ_2 = malloc(sizeR * sizeof(*PCQ_2)) ;
    VSPINT sep2Pcost = cost_p ;
    VSPINT sep2Qcost = cost_q ;
    VSPINT sep2Ccost = 0 ;

    for(i = 1 ; i <= sizeR ; i++)
    {
        v     = i-1 ;
        v_in  = 2*i - 1 ;
        v_out = v_in + 1 ;

        if( cut2[v_in] == 0 && cut2[v_out] == 0) // v goes in P
        {
            PCQ_2[v]   = 0 ;
            sep2Pcost += c[ Rdense[v] ] ;
        }
        else if(cut2[v_in] == 1 && cut2[v_out] == 1) // v goes in Q
        {
            PCQ_2[v]   = 1 ;
            sep2Qcost += c[ Rdense[v] ] ;
        }
        else // v goes in C
        {
            PCQ_2[v]   = 2 ;
            sep2Ccost += c[ Rdense[v] ] ;
        }
    }
    VSPINT sep2imb = ABS_VAL(sep2Pcost - sep2Qcost) ;

    // ========================================================================
    // Apply the best new separator
    // ========================================================================
    VSPINT original_sep_cost =   problem->graph->total_vtx_cost
                               - X->weight - Y->weight ;
    VSPINT original_imb = ABS_VAL(X->weight - Y->weight) ;

    // A partition is feasible if both shores satisfy the upper bound u
    int sep1feasible = (sep1Pcost <= problem->u) && (sep1Qcost <= problem->u) ;
    int sep2feasible = (sep2Pcost <= problem->u) && (sep2Qcost <= problem->u) ;

    // Take a feasible partition, if one exists. Both have the same separator
    // cost (the value of the max flow), so if both are feasible, take the
    // one with the smaller imbalance (partition 1 on ties)
    int    *best_PCQ   = NULL ;
    VSPINT  best_Ccost = 0 ;
    VSPINT  best_imb   = 0 ;
    if(sep1feasible && (!sep2feasible || sep1imb <= sep2imb))
    {
        best_PCQ   = PCQ_1 ;
        best_Ccost = sep1Ccost ;
        best_imb   = sep1imb ;
    }
    else if(sep2feasible)
    {
        best_PCQ   = PCQ_2 ;
        best_Ccost = sep2Ccost ;
        best_imb   = sep2imb ;
    }

    // Determine if the separator improved. The partition taken replaces
    // X, S, Y only if C(S') < C(S), or C(S') = C(S) and its imbalance is
    // smaller
    int moved = 0 ;
    if(best_PCQ != NULL
       && (    best_Ccost <  original_sep_cost
           || (best_Ccost == original_sep_cost && best_imb < original_imb)))
    {
        // Move each vertex w of R to its new side. P joins X (side 0), Q
        // joins Y (side 1), C goes to the separator (0 in both x and y).
        // Vertices outside of R keep their sides. For each of x and y whose
        // entry for w changes, update the entry, Az (add or subtract column
        // w of A) and the weight. Like Cover_Refine, this assumes that x and
        // y are binary (no fractional entry), so only int_part needs updating
        VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;
        State  *states[2] = {X, Y} ;
        VSPINT *z, w, want ;
        for(v = 0 ; v < sizeR ; v++)
        {
            w = Rdense[v] ;
            for(k = 0 ; k < 2 ; k++)
            {
                z    = states[k]->point->int_part ;
                want = (best_PCQ[v] == k) ;
                if(z[w] == want) continue ;

                moved = 1 ;
                z[w] = want ;
                Az_pm_col_of_A
                (states[k]->Az->int_part, w, want ? PLUS : MINUS, problem) ;
                states[k]->weight += want ? c[w] : -c[w] ;
            }
        }
    }

    // Free memory
    free(PCQ_1) ;
    free(PCQ_2) ;
    free(cut2) ;
    free(RTp) ;
    free(RTi) ;
    free(RTx) ;
    free(queue) ;
    free(cut1) ;
    free(ptrs) ;
    free(Rp) ;
    free(Ri) ;
    free(Rx) ;

    return moved ;
}

// Refine the partition represented by X and Y
void Max_Flow_Refine_2_cuts
(
    State   *X,
    State   *Y,
    Problem *problem
)
{
    igraph_vector_t flow ;
    VSPINT numA, vF, i, v, cost_p, cost_q ;
    int moved ;
    VSPINT  *Ai       = problem->graph->Ai ;
    VSPLONG *Ap       = problem->graph->Ap ;
    VSPINT   n        = problem->graph->n ;
    VSPINT  *sides[2] = {X->point->int_part, Y->point->int_part} ;
    VSPINT   wgts[3] ; // set at the start of every pass
    VSPINT   u        = problem->u ;
    VSPINT  *c        = problem->graph->c ;
    VSPINT   gamma    = problem->graph->max_vtx_cost ;
    
    VSPINT   deg_budget = 256 * sqrt(n) ;
    VSPINT  *Adense     = malloc(n * sizeof(*Adense)) ;
    VSPINT  *Asparse    = malloc(n * sizeof(*Asparse)) ;
    for(i = 0 ; i < n ; i++) Asparse[i] = EMPTY ;

    // F has at most 3 * numA + deg(A) arcs (see Generate_Max_Flow_Graph), and
    // numA <= n, deg(A) <= Ap[n]
    VSPINT  *Fi = malloc( (3 * (VSPLONG) n + Ap[n]) * sizeof(*Fi)) ;
    VSPINT  *Fx = malloc( (3 * (VSPLONG) n + Ap[n]) * sizeof(*Fx)) ;
    VSPLONG *Fp = malloc( (2 * n + 3) * sizeof(*Fp)) ;

    int *Rpart1  = malloc(n * sizeof(*Rpart1)) ;
    int *Rpart2  = malloc(n * sizeof(*Rpart2)) ;


    for(int num_its = 0 ; num_its < 3 ; num_its++)
    {
        // The weights change when vertices move, so refresh them every pass
        wgts[0] = X->weight ;
        wgts[1] = Y->weight ;
        wgts[2] = problem->graph->total_vtx_cost - X->weight - Y->weight ;

        // Skip the max flow if the separator is too small or too large for it
        // to pay off (see MF_MIN_SEP and MF_MAX_SEP_FACTOR at the top)
        if(wgts[2] < MF_MIN_SEP || wgts[2] > MF_MAX_SEP_FACTOR * sqrt(n)) break;

        // Get the region R. Quit if no region
        numA = Get_R(Ai, Ap, n, u, c, sides, wgts, deg_budget, Adense, Asparse);
        if(numA == 0) break ;


        // Get the flow graph corresponding to the area
        Generate_Max_Flow_Graph
        (Ai, Ap, n, c, gamma, Adense, Asparse, numA, sides, Fi, Fx, Fp) ;
        vF = 2*numA + 2 ; // number of vertices in s-t flow graph

        // Get the max flow
        flow = Get_Max_Flow(Fi, Fx, Fp, vF) ;

        // Cost of X \ A and Y \ A, which are contracted into p and q. This
        // must be done before Apply_Flow_Refinement moves the vertices of A
        cost_p = X->weight ;
        cost_q = Y->weight ;
        for(i = 0 ; i < numA ; i++)
        {
            v = Adense[i] ;
            if(sides[0][v]) cost_p -= c[v] ;
            if(sides[1][v]) cost_q -= c[v] ;
        }

        // Read off the two extreme min cuts, and apply the best one if it
        // improves the separator
        moved = Apply_Flow_Refinement
        (Fi, Fx, Fp, vF, gamma, flow, Adense, c, cost_p, cost_q, X, Y, problem) ;

        // Clean up for the next pass: free the flow, and reset Asparse to all
        // EMPTY, as Get_R assumes
        igraph_vector_destroy(&flow) ;
        for(i = 0 ; i < numA ; i++) Asparse[ Adense[i] ] = EMPTY ;

        // If nothing moved, X, S, Y are unchanged, so the next pass would
        // build the same A and find the same cuts. Stop
        if(!moved) break ;
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


// ============================================================================
// One cut version: igraph finds the min cut with the smallest Q itself, so we
// do not need the flow or the residual network
// ============================================================================

// Find a max flow on F with igraph, and set Qhat[x] = 1 for the nodes x of F
// from which q is reachable in the residual network (0 for all others). This
// is the minimum pq cut with the smallest Q-hat. igraph finds these nodes
// with a backward BFS from q (its "partition2" output). We do not ask igraph
// for the flow itself, which saves igraph a second phase that turns its
// preflow into a flow
static void Get_Smallest_Q_Cut
(
    // The Network
    VSPINT  *Fi,
    VSPINT  *Fx,
    VSPLONG *Fp,
    VSPINT   vF,   // number of vertices in F

    // Output: size vF
    int     *Qhat
)
{
    VSPINT  i ;
    VSPLONG p, end ;
    VSPINT  eF = Fp[vF] ; // number of edges in F

    // ========================================================================
    // Store the network edges and edge capacities in igraph vectors
    // ========================================================================
    igraph_vector_int_t edges ;
    igraph_vector_t capacities ;
    igraph_vector_int_init(&edges, 2*eF) ;
    igraph_vector_init(&capacities, eF) ;

    for(i = 0 ; i < vF ; i++)
    {
        end = Fp[i+1] ;
        for(p = Fp[i] ; p < end ; p++)
        {
            VECTOR(edges)[2*p]     = i ;
            VECTOR(edges)[2*p + 1] = Fi[p] ;
            VECTOR(capacities)[p]  = Fx[p] ;
        }
    }

    // ========================================================================
    // Create the igraph network and get the cut
    // ========================================================================
    igraph_t F ;
    igraph_create(&F, &edges, vF, IGRAPH_DIRECTED) ;

    igraph_vector_int_t partition2 ; // the nodes from which q is reachable
    igraph_vector_int_init(&partition2, 0) ;
    igraph_maxflow
    (&F, NULL, NULL, NULL, NULL, &partition2, 0, vF - 1, &capacities, NULL) ;

    for(i = 0 ; i < vF ; i++) Qhat[i] = 0 ;
    for(i = 0 ; i < igraph_vector_int_size(&partition2) ; i++)
    {
        Qhat[ VECTOR(partition2)[i] ] = 1 ;
    }

    // Free the igraph objects
    igraph_vector_int_destroy(&partition2) ;
    igraph_destroy(&F) ;
    igraph_vector_destroy(&capacities) ;
    igraph_vector_int_destroy(&edges) ;
}


// Read off the partition P, C, Q of R from the cut Qhat, and apply it if it
// is feasible and improves the separator. Same rule as in
// Apply_Flow_Refinement: the partition replaces X, S, Y only if C(S') < C(S),
// or C(S') = C(S) and its imbalance is smaller.
// Returns 1 if a vertex moved (X, S, Y changed), 0 if not
static int Apply_Smallest_Q_Cut
(
    int     *Qhat,   // Qhat[x] = 1 iff q is reachable from node x of F
    VSPINT   vF,     // number of vertices in F
    VSPINT  *Rdense, // the vertices of R
    VSPINT  *c,      // cost vector of the original graph
    VSPINT   cost_p, // cost of X \ R
    VSPINT   cost_q, // cost of Y \ R
    State   *X,
    State   *Y,
    Problem *problem
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, k, v, v_in, v_out, w, want ;
    VSPINT *z ;
    int     moved = 0 ;

    // ========================================================================
    // Figure out the cost of the shores given the separating partition P, C, Q
    // ========================================================================
    // Qhat marks the q side: v_in and v_out both reach q -> Q, neither
    // does -> P, only v_out does -> C (only v_in reaching q is impossible)
    //
    // In PCQ array, indexed by the position v of a vertex in R (Rdense):
    // 0 -- in P
    // 1 -- in Q
    // 2 -- in C
    VSPINT sizeR = (vF - 2) / 2 ;
    int   *PCQ   = malloc(sizeR * sizeof(*PCQ)) ;
    VSPINT Pcost = cost_p ;
    VSPINT Qcost = cost_q ;
    VSPINT Ccost = 0 ;

    for(i = 1 ; i <= sizeR ; i++)
    {
        v     = i-1 ;
        v_in  = 2*i - 1 ;
        v_out = v_in + 1 ;

        if( Qhat[v_in] == 0 && Qhat[v_out] == 0) // v goes in P
        {
            PCQ[v] = 0 ;
            Pcost += c[ Rdense[v] ] ;
        }
        else if(Qhat[v_in] == 1 && Qhat[v_out] == 1) // v goes in Q
        {
            PCQ[v] = 1 ;
            Qcost += c[ Rdense[v] ] ;
        }
        else // v goes in C
        {
            PCQ[v] = 2 ;
            Ccost += c[ Rdense[v] ] ;
        }
    }
    VSPINT imb = ABS_VAL(Pcost - Qcost) ;

    // ========================================================================
    // Apply the new separator if it is feasible and improves the old one
    // ========================================================================
    VSPINT original_sep_cost =   problem->graph->total_vtx_cost
                               - X->weight - Y->weight ;
    VSPINT original_imb = ABS_VAL(X->weight - Y->weight) ;
    int feasible = (Pcost <= problem->u) && (Qcost <= problem->u) ;

    if(feasible
       && (    Ccost <  original_sep_cost
           || (Ccost == original_sep_cost && imb < original_imb)))
    {
        // Move each vertex w of R to its new side. P joins X (side 0), Q
        // joins Y (side 1), C goes to the separator (0 in both x and y).
        // Vertices outside of R keep their sides. For each of x and y whose
        // entry for w changes, update the entry, Az (add or subtract column
        // w of A) and the weight. Like Cover_Refine, this assumes that x and
        // y are binary (no fractional entry), so only int_part needs updating
        VSP_ASSERT(X->point->ind_frac == EMPTY && Y->point->ind_frac == EMPTY) ;
        State *states[2] = {X, Y} ;
        for(v = 0 ; v < sizeR ; v++)
        {
            w = Rdense[v] ;
            for(k = 0 ; k < 2 ; k++)
            {
                z    = states[k]->point->int_part ;
                want = (PCQ[v] == k) ;
                if(z[w] == want) continue ;

                moved = 1 ;
                z[w] = want ;
                Az_pm_col_of_A
                (states[k]->Az->int_part, w, want ? PLUS : MINUS, problem) ;
                states[k]->weight += want ? c[w] : -c[w] ;
            }
        }
    }

    free(PCQ) ;
    return moved ;
}


// Refine the partition represented by X and Y, like Max_Flow_Refine_2_cuts,
// but using only the min cut with the smallest Q, which igraph returns
// directly. No flow and no residual network are needed
void Max_Flow_Refine_1_Cut
(
    State   *X,
    State   *Y,
    Problem *problem
)
{
    VSPINT numA, vF, i, v, cost_p, cost_q ;
    int moved ;
    VSPINT  *Ai       = problem->graph->Ai ;
    VSPLONG *Ap       = problem->graph->Ap ;
    VSPINT   n        = problem->graph->n ;
    VSPINT  *sides[2] = {X->point->int_part, Y->point->int_part} ;
    VSPINT   wgts[3] ; // set at the start of every pass
    VSPINT   u        = problem->u ;
    VSPINT  *c        = problem->graph->c ;
    VSPINT   gamma    = problem->graph->max_vtx_cost ;

    VSPINT   deg_budget = 256 * sqrt(n) ;
    VSPINT  *Adense     = malloc(n * sizeof(*Adense)) ;
    VSPINT  *Asparse    = malloc(n * sizeof(*Asparse)) ;
    for(i = 0 ; i < n ; i++) Asparse[i] = EMPTY ;

    // F has at most 3 * numA + deg(A) arcs (see Generate_Max_Flow_Graph), and
    // numA <= n, deg(A) <= Ap[n]
    VSPINT  *Fi = malloc( (3 * (VSPLONG) n + Ap[n]) * sizeof(*Fi)) ;
    VSPINT  *Fx = malloc( (3 * (VSPLONG) n + Ap[n]) * sizeof(*Fx)) ;
    VSPLONG *Fp = malloc( (2 * n + 3) * sizeof(*Fp)) ;

    // Qhat[x] = 1 iff q is reachable from node x of F. F has at most 2n + 2
    // nodes
    int *Qhat = malloc( (2 * n + 2) * sizeof(*Qhat)) ;


    for(int num_its = 0 ; num_its < 3 ; num_its++)
    {
        // The weights change when vertices move, so refresh them every pass
        wgts[0] = X->weight ;
        wgts[1] = Y->weight ;
        wgts[2] = problem->graph->total_vtx_cost - X->weight - Y->weight ;

        // Skip the max flow if the separator is too small or too large for it
        // to pay off (see MF_MIN_SEP and MF_MAX_SEP_FACTOR at the top)
        if(wgts[2] < MF_MIN_SEP || wgts[2] > MF_MAX_SEP_FACTOR * sqrt(n)) break;

        // Get the region R. Quit if no region
        numA = Get_R(Ai, Ap, n, u, c, sides, wgts, deg_budget, Adense, Asparse);
        if(numA == 0) break ;


        // Get the flow graph corresponding to the area
        Generate_Max_Flow_Graph
        (Ai, Ap, n, c, gamma, Adense, Asparse, numA, sides, Fi, Fx, Fp) ;
        vF = 2*numA + 2 ; // number of vertices in s-t flow graph

        // Get the min cut with the smallest Q
        Get_Smallest_Q_Cut(Fi, Fx, Fp, vF, Qhat) ;

        // Cost of X \ A and Y \ A, which are contracted into p and q. This
        // must be done before Apply_Smallest_Q_Cut moves the vertices of A
        cost_p = X->weight ;
        cost_q = Y->weight ;
        for(i = 0 ; i < numA ; i++)
        {
            v = Adense[i] ;
            if(sides[0][v]) cost_p -= c[v] ;
            if(sides[1][v]) cost_q -= c[v] ;
        }

        // Read off the partition from the cut, and apply it if it improves
        // the separator
        moved = Apply_Smallest_Q_Cut
        (Qhat, vF, Adense, c, cost_p, cost_q, X, Y, problem) ;

        // Clean up for the next pass: reset Asparse to all EMPTY, as Get_R
        // assumes
        for(i = 0 ; i < numA ; i++) Asparse[ Adense[i] ] = EMPTY ;

        // If nothing moved, X, S, Y are unchanged, so the next pass would
        // build the same A and find the same cut. Stop
        if(!moved) break ;
    }

    // Free memory
    free(Fi) ;
    free(Fx) ;
    free(Fp) ;
    free(Adense) ;
    free(Asparse) ;
    free(Qhat) ;
}
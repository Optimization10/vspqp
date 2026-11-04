#include <math.h>
#include <igraph.h>

// This file includes only igraph's headers, not VSPQP's, so that the two
// cannot clash. The prototype is repeated in solve_internal.h for
// band_refine.c.
void Band_Igraph_Max_Flow
(
    int nn, int ne, const int *to, long long *cap, int s, int t
) ;


// Maximum flow by igraph for Band_Refine. The band network
// (forward arcs e = 0, 2, 4, ..., reverse arc e+1, residual capacities in cap)
// is copied into an igraph graph, and the flow on each arc is written back
// into the residual capacities, so the cuts are read as after Edmonds-Karp.
// cap must hold the original capacities (nothing pushed yet).
void Band_Igraph_Max_Flow
(
    int        nn,
    int        ne,
    const int *to,
    long long *cap,
    int        s,
    int        t
)
{
    igraph_int_t        m = ne / 2, k ;
    igraph_vector_int_t edges ;
    igraph_vector_t     capv, flow ;
    igraph_t            g ;
    igraph_real_t       value ;

    igraph_vector_int_init(&edges, 2 * m) ;
    igraph_vector_init(&capv, m) ;
    igraph_vector_init(&flow, 0) ;
    for(k = 0 ; k < m ; k++)
    {
        VECTOR(edges)[2*k]   = to[2*k+1] ;   // tail of arc 2k
        VECTOR(edges)[2*k+1] = to[2*k] ;     // head of arc 2k
        VECTOR(capv)[k]      = (double) cap[2*k] ;
    }

    igraph_create(&g, &edges, nn, IGRAPH_DIRECTED) ;
    igraph_maxflow(&g, &value, &flow, NULL, NULL, NULL, s, t, &capv, NULL) ;

    for(k = 0 ; k < m ; k++)
    {
        long long f = (long long) llround(VECTOR(flow)[k]) ;
        cap[2*k]   -= f ;
        cap[2*k+1] += f ;
    }

    igraph_destroy(&g) ;
    igraph_vector_destroy(&flow) ;
    igraph_vector_destroy(&capv) ;
    igraph_vector_int_destroy(&edges) ;
}

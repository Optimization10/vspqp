#include "shared_api.h"
#include "coarsen_internal.h"

/* Get all matches of vertex. The "matches" array will include vertex
   itself at index 0 */
VSPINT* VSP_get_all_matches (
    VSPINT    vertex,
    VSPINT *matching
) {
    static long long unsigned count = 0 ;
    count += 1 ;
    VSPINT i ;
    VSPINT* matches = malloc (4 * sizeof(VSPINT)) ;
    for (i = 0 ; i < 4 ; i++) matches [i] = EMPTY ;

    matches [0] = vertex ;
    matches [1] = VSP_get_match(matches [0] , matching) ;

    if (matches [1] == vertex) {
        matches [1] = EMPTY ;
    }
    else {
        matches [2] = VSP_get_match(matches [1], matching) ;
        if (matches [2] == vertex) {
            matches [2] = EMPTY ;
        }
        else {
            if ( VSP_get_match(matches [2], matching) != vertex ) {
                VSP_error("Error with matching", DBG_INFO) ;
            }

        }
    }

    return matches ;
}


/* Get the edge weight between vertex_a and vertex_b */
VSPINT VSP_get_edge_weight (
    VSPINT      vertex_a,
    VSPINT      vertex_b,
    VSPINT           *Ai,
    VSPLONG *Ap,
    VSPINT *Ae
) {

    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    VSPLONG p, end ;
    VSPINT i, weight ;

    end = Ap [vertex_a + 1] ;

    weight = 0 ;
    for ( p = Ap [vertex_a] ; p < end ; p++ ) {
        i = Ai [p] ;
        if ( i == vertex_b ) {
            weight = Ae [p] ;
            break ;
        } 
    }

    #ifdef VSP_DEBUG

        /* Get the weight of the "reverse" edge */
        VSPINT reverse_weight = FALSE ;

        end = Ap [vertex_b + 1] ;

        for ( p = Ap [vertex_b] ; p < end ; p++ ) {
            i = Ai [p] ;
            if ( i == vertex_a ) {
                reverse_weight = Ae [p] ;
                break ;
            }
        }
        if (weight != reverse_weight) {
            VSP_error("Edge weights not symmetric", DBG_INFO) ;
        }

    #endif

    return weight ;
}


/* Given two arrays of vertices, determine if there exists an edge connecting
   any vertex from the first array to any vertex of the second array          
   
   It is assumed that both of the arrays have the form
   { v1, v2, ..., vk, EMPTY, ..., EMPTY }
   It is okay if one, or both, are completely empty, or completely full       */

VSPINT VSP_check_for_edge_between_arrays_of_vertices (
    VSPINT *arr1,
    VSPINT *arr2,
    VSPINT  size,
    VSPINT   *Ai,
    VSPLONG *Ap
) {
    VSPINT r, s, v1, v2 ;
    VSPINT there_is_an_edge = FALSE ;

    WARN_IF_NOT_DEBUG_MODE ;

    v1 = arr1 [ r = 0 ] ;

    while ( v1 != EMPTY && r < size ) {

        v2 = arr2 [ s = 0 ] ;

        while (v2 != EMPTY && s < size ) {

            if ( VSP_there_is_an_edge(v1, v2, Ai, Ap) ) {
                there_is_an_edge = TRUE ;
            }
            v2 = arr2 [ ++s ] ;
        }

        v1 = arr1 [ ++ r ] ;
    }

    return there_is_an_edge ;

}

/* Check if there is an edge between vertex_a and vertex_b */
VSPINT VSP_there_is_an_edge (
    VSPINT vertex_a,
    VSPINT vertex_b,
    VSPINT      *Ai,
    VSPLONG *Ap
) {

    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    VSPLONG p, end ;
    VSPINT i, there_is_an_edge = FALSE ;

    end = Ap [vertex_a + 1] ;

    /* Check that vertex_b is among the neighbors of vertex_a */
    for ( p = Ap [vertex_a] ; p < end ; p++ ) {
        i = Ai [p] ;
        if ( i == vertex_b ) {
            there_is_an_edge = TRUE ;
            break ;
        } 
    }

    #ifdef VSP_DEBUG

        /* Check that vertex_a is among the neighbors of vertex_b */
        VSPINT reverse_edge = FALSE ;

        end = Ap [vertex_b + 1] ;

        for ( p = Ap [vertex_b] ; p < end ; p++ ) {
            i = Ai [p] ;
            if ( i == vertex_a ) {
                reverse_edge = TRUE ;
                break ;
            }
        }
        if (reverse_edge != there_is_an_edge) {
            VSP_error("Matrix not symmetric", DBG_INFO) ;
        }

    #endif

    return there_is_an_edge ;
}

/* Given a problem and a vertex in the coarse graph */
VSPINT* VSP_get_preimage_of_vertex
(
    VSPINT vertex_coarse,         /* A vertex in the coarse problem */
    VSPINT     *matching,
    VSPINT  *invmatchmap
)
{
    VSPINT i ;

    VSPINT* preimage = malloc ( 4 * sizeof(VSPINT) ) ;
    for (i = 0 ; i < 4 ; i++) preimage [i] = EMPTY ;

    preimage [0] = invmatchmap [vertex_coarse] ;
    preimage [1] = VSP_get_match(preimage [0], matching);

    if ( preimage [1] == preimage [0] ) {
        preimage [1] = EMPTY ;
    }
    else {
        preimage [2] = VSP_get_match(preimage [1], matching) ;

        if (preimage [2] == preimage [0] ) {
            preimage [2] = EMPTY ;
        }
    }

    return preimage ;
}

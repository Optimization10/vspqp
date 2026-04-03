#include "core.h"
#include "shared_api.h"
#include "coarsen_internal.h"
#include "shared_macros.h"
#include "shared_types.h"

/* Matching function declarations */
void Heavy_Edge_Matching(VSP_problem* problem, VSP_options* options);
void VSP_matching_SR_HEM(VSP_problem* problem, VSP_options* options);




void Match_Check
(
    VSPINT                     n,
    VSPINT                    cn,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    VSP_match_type    *matchtype
)
{
    static long long unsigned count = 0 ;
    count ++ ;

    VSPINT i, k0, k1, k2, k3 ;
    VSPINT matched_to_invalid_vertex = FALSE ;
    VSPINT wrong_matchtype = FALSE ;


    /* Check that matchmap composed with invmatchmap is the identiy */
    for (i = 0 ; i < cn ; i++) {
        if ( matchmap [ invmatchmap [i] ] != i ) {
            VSP_error("Error in matching.", DBG_INFO) ;
        }
    }


    /* Check that each vertex in the fine graph is matched with something */
    for (i = 0 ; i < n ; i++)
    {
        if ( matching [i] <= 0 || matching [i] >= n + 1 )
        {
            matched_to_invalid_vertex = TRUE ;
            break ;
        }
        if
        (
            // Matchtypes for Tim Davis coarsening
            matchtype [i] != matchtype_standard  &&
            matchtype [i] != matchtype_brotherly &&
            matchtype [i] != matchtype_community &&
            matchtype [i] != matchtype_orphan    &&

            // Matchtypes for low_gamma_matching
            matchtype [i] != matchtype_HEM       &&
            matchtype [i] != matchtype_singleton &&
            matchtype [i] != matchtype_selfMatch &&
            matchtype [i] != matchtype_3way       
        )
        {
            wrong_matchtype = TRUE ;
            break ;
        }
    }
    if (matched_to_invalid_vertex)
    {
        VSP_error("Matched to invalid vertex", DBG_INFO) ;
    }
    
    if (wrong_matchtype)
    {
        VSP_error("Invalid matchtype", DBG_INFO) ;
    }


    /* Check that each member of each linked list is mapped to the same
        vertex in the coarse Graph-> */
    for( k0 = 0 ; k0 < n ; k0++ )
    {
        wrong_matchtype = FALSE ;

        k1 = matching [k0] - 1 ; /* k0 points to k1 */
        k2 = matching [k1] - 1 ; /* k1 points to k2 */
        k3 = matching [k2] - 1 ; /* k2 points to k3 */

        /* If k0 is matched with itself */
        if ( k1 == k0 )
        {
            if (matchtype [k0] != matchtype_orphan && matchtype [k0] != matchtype_selfMatch && matchtype [k0] != matchtype_singleton)
            {
                VSP_error("Matchtype incorrect", DBG_INFO) ;
            }
        }

        /* If k0 is matched with precisely 1 other vertex */
        else if ( k2 == k0 )
        {
            if ( matchmap [k0] != matchmap [k1] )
            {
                VSP_error("Error in matchmap", DBG_INFO) ;
            }
            if ( matchtype [k0] == matchtype_orphan ||
                    matchtype [k1] == matchtype_orphan    )
            {
                VSP_error("Error with matchtype", DBG_INFO) ;
            }
        }
        /* If k0 is matched with precisely 2 other vertices */
        else if ( k3 == k0 )
        {
            if ( matchmap [k0] != matchmap [k1] || 
                    matchmap [k0] != matchmap [k2]     )
            {
                VSP_error("Error in matchmap", DBG_INFO) ;
            }
            if ( matchtype [k0] == matchtype_orphan ||
                    matchtype [k1] == matchtype_orphan ||
                    matchtype [k2] == matchtype_orphan    )
            {
                VSP_error("Error with matchtype", DBG_INFO) ;
            }
        }
        else {
            /* Code should never reach here */
            VSP_error("Error with matching ll", DBG_INFO) ;
        }
    }
}



VSPINT VSP_get_match(VSPINT vertex, VSPINT* matching)
{
    return matching[vertex] - 1;
}


bool Is_Matched(VSPINT vertex, VSPINT* matching)
{
    return matching[vertex] ;
}

bool IsMatchedMatchtype
(
    VSPINT vertex,
    VSPINT* matchingArray,
    VSP_match_type *matchtypeArray,
    VSP_match_type desiredType
)
{
    return Is_Matched(vertex, matchingArray) && (matchtypeArray[vertex] == desiredType) ;
}
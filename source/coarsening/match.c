#include "shared_api.h"
#include "coarsen_internal.h"
#include "shared_macros.h"
#include "shared_types.h"

/* Matching function declarations */
void Heavy_Edge_Matching(VSP_problem* problem, VSP_options* options);
void VSP_matching_SR_HEM(VSP_problem* problem, VSP_options* options);


void VSP_create_match
(
    VSP_problem* problem, 
    VSPINT vertexA, VSPINT vertexB, 
    VSP_match_type matchtype
) ;

void VSP_matching_cleanup(VSP_problem* problem, VSP_options* options);
void VSP_create_community_match
(
    VSP_problem* problem, 
    VSPINT vertexA, 
    VSPINT vertexB, 
    VSP_match_type matchtype
) ;


void VSP_match(VSP_problem* problem, VSP_options* options)
{
    #ifdef VSP_DEBUG_COARSENING
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    switch (options->matching_strategy)
    {
        case HEM:
            Heavy_Edge_Matching(problem, options);
            VSP_matching_cleanup(problem, options) ;
            break;
        case HEMSR:
            Heavy_Edge_Matching(problem, options);
            VSP_matching_SR_HEM(problem, options);
            VSP_matching_cleanup(problem, options) ;
            break;
        // case LOW_GAMMA_MATCHING:
            // Low_Gamma_Matching(problem, options, 0.7) ;
            break ;
        default:
            printf("Must provide a matching strategy\n");
    }

    #ifdef VSP_DEBUG_COARSENING
        Match_Check
        (
            problem->Graph->n,
            problem->cn,
            problem->matching,
            problem->matchmap,
            problem->invmatchmap,
            problem->matchtype
        ) ;
    #endif

}


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


/* ========================================================================== */
/* === HEM Functions ======================================================== */
/* ========================================================================== */

void Heavy_Edge_Matching(VSP_problem* problem, VSP_options* options)
{
    /* Initialize local variables */
    VSPINT k = 0, p = 0, weight = 0;

    VSPINT ismatched = -1;
    VSPINT heaviestneighbor = -1;
    VSPINT heaviestweight = -1;


    /* Graph data */
    VSPINT n   = problem->Graph->n;
    VSPINT* Ai = problem->Graph->Ai;
    VSPLONG* Ap = problem->Graph->Ap;
    const VSPINT isoriginal = Is_Finest_Problem(problem);

    /* Matching data*/
    VSPINT* edgeweights = problem->Graph->Ae;
    VSPINT* matching = problem->matching;

    for (k = 0; k < n; k++)
    {
        /* Note: the +1 we do in matching is so that the following line works*/
        ismatched = matching[k];
        if (ismatched) continue ;

        /* Reset heaviest neighbor and its weight */
        heaviestneighbor = -1;
        heaviestweight = -1;

        for (p = Ap[k]; p < Ap[k+1]; p++)
        {
            VSPINT neighbor = Ai[p];
            /* Consider only unmatched neighbors */
            ismatched = matching[neighbor];
            if (ismatched) continue ;

            /* Keep track of heaviest */
            weight = isoriginal ? 1 : edgeweights[p];
            if (weight > heaviestweight)
            {
                heaviestweight   = weight;
                heaviestneighbor = neighbor;
            }
        }

        /* Match to the heaviest */
        if (heaviestneighbor != -1)
        {
            /* Match k and heaviest neighbor */
            VSP_create_match(problem, k, heaviestneighbor, matchtype_standard);
        }
    }
}

/* Pretty sure that HEM + this takes care of everything except for singletons */
void VSP_matching_SR_HEM(VSP_problem* problem, VSP_options* options)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    VSPINT n = problem->Graph->n;
    VSPLONG* Gp = problem->Graph->Ap;
    VSPINT* Gi = problem->Graph->Ai;
    VSPINT* Gx = problem->Graph->Ae;
    const bool isFinestProb = Is_Finest_Problem(problem) ;

    for(VSPINT k = 0; k < n; k++)
    {
        /* Consider only unmatched vertices */
        if (problem->matching[k])
            {continue;}
        VSPINT heaviestneighbor = -1;
        VSPINT heaviestweight = -1;

        /* Go through all the neighbors of k */
        for (VSPLONG p = Gp[k]; p < Gp[k+1]; p++)
        {
            VSPINT neighbor = Gi[p];
            #ifdef VSP_DEBUG
                if(problem->matching[neighbor] == 0)
                {
                    VSP_error("Neighbor is unmatched", DBG_INFO);
                }
            #endif

            /* Keep track of heaviest */
            VSPINT x = isFinestProb ? 1 : Gx[p];
            if (x > heaviestweight)
            {
                heaviestweight = x;
                heaviestneighbor = neighbor;
            }
        }

        /* If we found a heaviest neighbor, begin resolving matches */
        /* Here is where we assume that regular HEM ran already. If it didn't,
           we would match k with the heaviest neighbor. Instead, we're assuming
           that heaviest neighbor is matched already, and we're doing brotherly
           matching through it */
        /* Do brotherly matching with heaviestneighbor as the center */
        if (heaviestneighbor != -1)
        {
            /* v, neighbor are both unmatched indices in N(heaviestneighbor)
               We lost k, but k will get resolved since it is in
               N(heaviestneighbor) */
            VSPINT v = -1; 
            for (VSPLONG p = Gp[heaviestneighbor]; p < Gp[heaviestneighbor + 1];
                                                                        p++)
            {
                VSPINT neighbor = Gi[p];
                if (problem->matching[neighbor]) continue ;

                if (v == -1) v = neighbor ;
                else
                {
                    VSP_create_match(problem, v, neighbor, matchtype_brotherly);
                    v = -1;
                }
            }

            /* If we had a vertex left over */
            if (v != -1)
            {
                if(options->do_community_match)
                {
                    VSP_create_community_match
                    (
                        problem,
                        heaviestneighbor,
                        v,
                        matchtype_community
                    );
                }
                else
                {
                    VSP_create_match(problem, v, v, matchtype_orphan);
                }
            }
        }
    }

    #ifdef VSP_DEBUG
        /* If community match is enabled, every vertex should be matched
           at this point, except singletons */
        if(options->do_community_match)
        {
            for(VSPINT i = 0; i < n; i++)
            {
                if(problem->matching[i] == 0 && problem->Graph->Ap[i+1] - problem->Graph->Ap[i] > 0)
                {
                    VSP_error("Vertex is unmatched", DBG_INFO);
                    break ;
                }
            }
        }
    #endif
}

void VSP_create_match
(
    VSP_problem*     problem, 
    VSPINT           vertexA,
    VSPINT           vertexB, 
    VSP_match_type matchtype
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        if(Is_Matched(vertexA, problem->matching) || Is_Matched(vertexB, problem->matching))
        {
            ERR("Trying to match an already matched vertex") ;
        }
    #endif
    VSPINT cn = problem->cn;
    
    /* Create the match */
    problem->matching[vertexA] = vertexB + 1;
    problem->matching[vertexB] = vertexA + 1;
    problem->invmatchmap[cn] = vertexA;
    problem->matchtype[vertexA] = matchtype;
    problem->matchtype[vertexB] = matchtype;
    problem->matchmap[vertexA] = cn;
    problem->matchmap[vertexB] = cn;
    problem->cn += 1;
}


/* Match singletons in pairs. If there is an odd number of singletons, match
   the last singleton to itself, or do a community match if that is enabled.

   For unmatched non-singletons, match them to themselves. If community matching
   is enabled, attempt to match them to a neighbor that is not community
   matched. If no such neighbor exists, then community match them to a neighbor.
*/
void VSP_matching_cleanup(VSP_problem* problem, VSP_options* options)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT k, ismatched, degree, singleton, p;
    VSPINT n = problem->Graph->n;
    VSPINT* Ap = problem->Graph->Ap;
    VSPINT* matching = problem->matching;

    /* Match unmatched vertices to themselves */
    for (k = 0; k < n; k++)
    {
        /* Match unmatched vertices to themselves */
        ismatched = matching[k];
        if ( !(ismatched) )
        {
            degree = Ap[k+1] - Ap[k];
            /* If k is a singleton */
            if (degree == 0)
            {
                if (problem->singleton == -1)
                {
                    problem->singleton = k;
                }
                else
                {
                    singleton = problem->singleton;
                    VSP_create_match(problem, k, singleton, matchtype_standard);
                    problem->singleton = -1;
                }
            }
            /* Else, if k is not a singleton */
            else
            {
                if(options->do_community_match)
                {
                    #ifdef VSP_DEBUG
                        /* In my code, this match cleanup function comes after
                        brotherly matching. If brotherly matching happened
                        already WITH community matching enabled, then there
                        shouldn't be any unmatched non-singletons */
                        VSP_error("Vertex should have been matched", DBG_INFO);
                    #endif

                    /* If possible, match k with a neighbor which doesn't have
                       matchtype "community matching"                         */
                    for (p = Ap [k] ; p < Ap [k+1] ; p++ )
                    {
                        VSPINT node = problem->Graph->Ai[p] ;
                        if(problem->matchtype[ node ] != matchtype_community)
                        {
                            break ;
                        }
                    }
                    VSP_create_community_match
                    (
                        problem, 
                        problem->Graph->Ai[p], 
                        k, 
                        matchtype_community
                    );
                }
                else
                {
                    VSP_create_match(problem, k, k, matchtype_orphan) ;
                }
            }
        }
    }

    /* If we have a leftover singleton */
    if(problem->singleton != -1)
    {
        VSPINT k = problem->singleton;
        if(options->do_community_match)
        {
            /* Attempt to community match the singleton with some vertex that
               does not have matchtype community */
            VSPINT i;
            for(i = 0; i < problem->Graph->n; i++)
            {
                /* Why i != k? If 0 is the one and only singleton, that would
                   be a serious bug which leads to invalid memory access.
                   matching[0] would be garbage, so we may choose to match 0
                   with 0. But then in create_community_match, we get the match
                   of 0, assuming that 0 is matched. This results in -1. Then,
                   we get the match of -1, which is where the real problem is.*/
                if(problem->matchtype[i] != matchtype_community && i != k)
                {
                    break;
                }

            }
            VSP_create_community_match(problem, i, k, matchtype_community) ;
        }
        else
        {
            VSP_create_match(problem, k, k, matchtype_orphan);
        }
    }
}

void VSP_create_community_match
(
    VSP_problem* problem, 
    VSPINT vertexA,          /* Assumes that vertexA is matched with SOMETHING*/
    VSPINT vertexB, 
    VSP_match_type matchtype
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    VSPINT* matching = problem->matching ;

    /* Get all the matches of vertex A */
    VSPINT vm[4] = {-1, -1, -1, -1};
    vm[0] = vertexA;
    vm[1] = VSP_get_match(vm[0], matching);
    #ifdef VSP_DEBUG
        if(vm[1] < 0 || vm[1] > problem->Graph->n - 1)
        {
            VSP_error("Invalid access", DBG_INFO) ;
        }
    #endif
    vm[2] = VSP_get_match(vm[1], matching);
    #ifdef VSP_DEBUG
        if(vm[2] < 0 || vm[2] > problem->Graph->n - 1)
        {
            VSP_error("Invalid access", DBG_INFO) ;
        }
    #endif
    vm[3] = VSP_get_match(vm[2], matching);
    #ifdef VSP_DEBUG
        if(vm[3] < 0 || vm[3] > problem->Graph->n - 1)
        {
            VSP_error("Invalid access", DBG_INFO) ;
        }
    #endif

    /* Determine if vertex A is part of a 3-way match */
    VSPINT is3way = (vm[0] == vm[3]);
    
    if (is3way)
    {
        /* Here, we steal one vertex from the 3way matching. But it is possible
           that invmatchmap pointed to this vertex we stole. If it did, we
           must change the invmatchmap */
        
        if ( problem->invmatchmap[ problem->matchmap[ vm[0] ] ] == vm[2] ) {
            problem->invmatchmap[ problem->matchmap[ vm[0] ] ] = vm [0] ;
        }

        problem->matching[vm[1]] = vertexA + 1;
        VSP_create_match(problem, vm[2], vertexB, matchtype);
    }
    else
    {
        problem->matching[vertexB] = problem->matching[vertexA];
        problem->matching[vertexA] = vertexB + 1;
        problem->matchmap[vertexB] = problem->matchmap[vertexA];
        problem->matchtype[vertexB] = matchtype;
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
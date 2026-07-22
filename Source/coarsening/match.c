#include "core.h"
#include "coarsen_internal.h"
#include "data_structures.h"
#include "utils.h"


static void Match_Check
(
    VSPINT    n, // fine n
    Matching *M
)
{
    static long long unsigned count = 0 ;
    count ++ ;

    VSPINT cn = M->cn ;
    VSPINT *matching = M->match ;
    VSPINT *matchmap = M->matchmap ;
    VSPINT *invmatchmap = M->invmatchmap ;
    match_type *matchtype = M->matchtype ;

    VSPINT i, k0, k1, k2, k3 ;
    VSPINT matched_to_invalid_vertex = FALSE ;
    VSPINT wrong_matchtype = FALSE ;


    /* Check that matchmap composed with invmatchmap is the identiy */
    for (i = 0 ; i < cn ; i++) {
        if ( matchmap [ invmatchmap [i] ] != i ) {
            Error("Error in matching.", DBG_INFO) ;
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
        Error("Matched to invalid vertex", DBG_INFO) ;
    }
    
    if (wrong_matchtype)
    {
        Error("Invalid matchtype", DBG_INFO) ;
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
                Error("Matchtype incorrect", DBG_INFO) ;
            }
        }

        /* If k0 is matched with precisely 1 other vertex */
        else if ( k2 == k0 )
        {
            if ( matchmap [k0] != matchmap [k1] )
            {
                Error("Error in matchmap", DBG_INFO) ;
            }
            if ( matchtype [k0] == matchtype_orphan ||
                 matchtype [k1] == matchtype_orphan     )
            {
                Error("Error with matchtype", DBG_INFO) ;
            }
        }
        /* If k0 is matched with precisely 2 other vertices */
        else if ( k3 == k0 )
        {
            if ( matchmap [k0] != matchmap [k1] || 
                    matchmap [k0] != matchmap [k2]     )
            {
                Error("Error in matchmap", DBG_INFO) ;
            }
            if ( matchtype [k0] == matchtype_orphan ||
                    matchtype [k1] == matchtype_orphan ||
                    matchtype [k2] == matchtype_orphan    )
            {
                Error("Error with matchtype", DBG_INFO) ;
            }
        }
        else {
            /* Code should never reach here */
            Error("Error with matching ll", DBG_INFO) ;
        }
    }
}



VSPINT get_match(VSPINT vertex, VSPINT* match)
{
    return match[vertex] - 1;
}


static bool Is_Matched(VSPINT vertex, VSPINT* matching)
{
    return matching[vertex] ;
}

static bool IsMatchedMatchtype
(
    VSPINT vertex,
    VSPINT* matchingArray,
    match_type *matchtypeArray,
    match_type desiredType
)
{
    return Is_Matched(vertex, matchingArray) && (matchtypeArray[vertex] == desiredType) ;
}


static void Additional_Match_Check
(
    Graph             *Graph,
    VSPINT                   lvl,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    match_type    *matchtype,
    VSPINT                    cn
) ;

static void Create_Match_2way
(
    VSPINT                       i,
    VSPINT                       j,
    VSPINT*               matching,
    VSPINT*            invmatchmap,
    match_type*      matchtype_arr,
    VSPINT*               matchmap,
    VSPINT*                 cn_ptr,
    match_type           mtype
)
{
    SET_FUNCTION_COUNTER;
    #ifdef DEBUG
        if(Is_Matched(i, matching) || Is_Matched(j, matching))
        {
            ERR("Trying to match an already matched vertex");
        }
    #endif
    VSPINT cn = *cn_ptr;
    
    /* Create the match */
    matching[i] = j + 1;
    matching[j] = i + 1;
    invmatchmap[cn] = i;
    matchtype_arr[i] = mtype;
    matchtype_arr[j] = mtype;
    matchmap[i] = cn;
    matchmap[j] = cn;
    *cn_ptr += 1;
}

static void Create_Match_3way_HEM
(
    VSPINT             i, // HEM
    VSPINT             j, // HEM
    VSPINT             k, // Unmatched
    VSPINT*            matching,
    VSPINT*            matchmap,
    VSPINT*            invmatchmap,
    match_type*    matchtype_arr
)
{

    SET_FUNCTION_COUNTER;
    #ifdef DEBUG
        if(Is_Matched(k, matching))
        {
            ERR("Trying to match an already matched vertex");
        }
        if(!IsMatchedMatchtype(i, matching, matchtype_arr, matchtype_HEM))
        {
            ERR("Non-HEM node in 3-way match");
        }
        if(!IsMatchedMatchtype(j, matching, matchtype_arr, matchtype_HEM))
        {
            ERR("Non-HEM node in 3-way match");
        }
    #endif


    matching[i] = j + 1;
    matching[j] = k + 1;
    matching[k] = i + 1;

    matchtype_arr[i] = matchtype_3way;
    matchtype_arr[j] = matchtype_3way;
    matchtype_arr[k] = matchtype_3way;

    VSPINT cv = matchmap[i] ;
    VSP_ASSERT(cv == matchmap[j]) ;
    matchmap[k] = cv ;

    VSP_ASSERT(invmatchmap[cv] == i || invmatchmap[cv] == j) ;
    
    return;
}

static void Create_Match_3way_Singleton
(
    VSPINT             i, // singleton, unmatched
    VSPINT             j, // singleton, unmatched
    VSPINT             k, // singleton, unmatched
    VSPINT*            matching,
    VSPINT*            matchmap,
    VSPINT*            invmatchmap,
    match_type*    matchtype_arr,
    VSPINT*            cn_ptr
)
{

    SET_FUNCTION_COUNTER;
    #ifdef DEBUG
        if(Is_Matched(i, matching))
        {
            ERR("Trying to match an already matched vertex");
        }

        if(Is_Matched(j, matching))
        {
            ERR("Trying to match an already matched vertex");
        }

        if(Is_Matched(k, matching))
        {
            ERR("Trying to match an already matched vertex");
        }

    #endif


    matching[i] = j + 1;
    matching[j] = k + 1;
    matching[k] = i + 1;

    matchtype_arr[i] = matchtype_singleton ;
    matchtype_arr[j] = matchtype_singleton ;
    matchtype_arr[k] = matchtype_singleton ;

    matchmap[i] = *cn_ptr ;
    matchmap[j] = *cn_ptr ;
    matchmap[k] = *cn_ptr ;
    
    invmatchmap[*cn_ptr] = i ;

    *cn_ptr += 1 ;
    return;
}


// ASSUME: old is in a 2-way brotherly match with some vertex i
// DO: Unmatch old and brotherly match new with i
static void Steal_Match
(
    VSPINT                     old,
    VSPINT                     new,
    VSPINT*               matching,
    VSPINT*            invmatchmap,
    match_type*      matchtype,
    VSPINT*               matchmap
)
{
    VSPINT i = get_match(old, matching) ;
 
    // Update matching
    matching[old] = 0 ;
    matchtype[old] = matchtype_unmatched ;
    matching[new] = i + 1 ;
    matching[i]   = new + 1 ;

    // Update invmatchmap
    VSPINT cv = matchmap[i] ;
    invmatchmap[cv] = i ;

    // Update matchmap
    matchmap[new] = cv ; 

    // Update matchtype
    matchtype[new] = matchtype_brotherly ;
}


// {i, j} = a HEM pair
// neighs = the adjacency list of i concatenated with the adjacency list of j
// numNeighs = size of neighs

// The following are invariants for the two loops visiting neighs:
// p = an index into neighs
// perfect  = the only unmatched entry of neighs[0...p-1] with cost 2^lvl. Or EMPTY
// light    = the only unmatched entry of neighs[0...p-1] with cost < 2^lvl. Or EMPTY
// lightest = the lowest-cost entry in neighs[0...p-1] that was unmatched before
//            starting the loop

static bool Is_Perfect_Vtx(VSPINT i, VSPINT *c, VSPINT pow2l)
{
    return (c[i] == pow2l) ;
}




static void Match_Singletons
(
    VSPINT           *heap_sings,
    VSPINT          *lookup_heap,
    VSPINT         *n_heap_sings,
    VSPINT                    *c,
    VSPINT          gamma_coarse,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    match_type    *matchtype,
    VSPINT                   *cn,
    VSPINT            *n_matches,
    VSPINT          *n_unmatched
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, j, k, ck ;
    MinHeap_Heapify_v0(heap_sings, lookup_heap, *n_heap_sings, c) ;

    while(*n_heap_sings > 0)
    {
        MinHeap_Get_Top_Val_v0(heap_sings, lookup_heap,n_heap_sings, c,&i);
        MinHeap_Get_Top_Val_v0(heap_sings, lookup_heap,n_heap_sings, c,&j);
        MinHeap_See_Top_Key_v0(heap_sings, c, *n_heap_sings, &ck) ;

        VSP_ASSERT(i != EMPTY) ;
        if(j != EMPTY && ck != EMPTY && c[i] + c[j] + ck <= gamma_coarse)
        {
            // extract k
            MinHeap_Get_Top_Val_v0(heap_sings, lookup_heap,n_heap_sings, c,&k);
            Create_Match_3way_Singleton
            (
                i, j, k, matching, matchmap, invmatchmap, matchtype, cn
            ) ;
            (*n_matches)++ ;
            (*n_unmatched) -= 3 ;
        }
        else if(j != EMPTY && c[i] + c[j] <= gamma_coarse)
        {
            // Match i with j
            Create_Match_2way
            (
                i, j, matching, invmatchmap, matchtype,
                matchmap, cn, matchtype_singleton
            ) ;
            (*n_matches)++ ;
            (*n_unmatched) -= 2 ;
        }
        else
        {
            // Match i with itself
            Create_Match_2way
            (
                i, i, matching, invmatchmap, matchtype,
                matchmap, cn, matchtype_singleton
            ) ;
            (*n_matches)++ ;
            (*n_unmatched) -= 1 ;

            if(j != EMPTY)
            {
                // Match j with itself
                Create_Match_2way
                (
                    j, j, matching, invmatchmap, matchtype,
                    matchmap, cn, matchtype_singleton
                ) ;
                (*n_matches)++ ;
                (*n_unmatched) -= 1 ;
            }
        }
    }
}

// Can handle states where not all vertices are matched yet
static void Check_cn
(
    VSPINT         n,
    VSPINT        cn, 
    VSPINT *matching
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, j, k ;
    VSPINT *seen = calloc(n, sizeof(VSPINT)) ;
    VSPINT cn_local = 0 ;

    for(i = 0 ; i < n ; i++)
    {
        if(seen[i] || !Is_Matched(i, matching)) continue ;

        cn_local++ ;
        j = get_match(i, matching) ;
        k = get_match(j, matching) ;
        seen[i] = 1 ;
        seen[j] = 1 ;
        seen[k] = 1 ;
    }

    #ifdef DEBUG
        if(cn != cn_local)
        {
            ERR("Wrong number of coarse vertices") ;
        }
    #endif

    free(seen) ;
}


// Check the following:
// A vertex and its match(es) all have the same matchtype

// If i has matchtype HEM, then it's matched to precisely one other vertex,
// and this vertex is a neihbor

// If i has matchtype brotherly, then it is matched to precisely one other
// vertex and this vertex is NOT a neighbor

// If i is in a 3-way match, then it is matched to precisely 2 other vertices,
// and it neighbors at least one of them

// cn is correct

static void Additional_Match_Check
(
    Graph                 *Graph,
    VSPINT                   lvl,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    match_type        *matchtype,
    VSPINT                    cn
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT i, j, k ;
    VSPLONG p, end ;

    // Graph info
    VSPINT   n   = Graph->n ;
    VSPLONG  nz  = Graph->nz ;
    VSPINT  *Ai  = Graph->Ai ;
    VSPINT  *Ae  = Graph->Ae ;
    VSPLONG *Ap  = Graph->Ap ;
    VSPINT  *c   = Graph->c ;


    // If i has matchtype HEM, then it's matched to precisely one other vertex,
    // whose matchtype is also HEM, and these vertices are neighbors
    for(i = 0 ; i < n ; i++)
    {
        j = get_match(i, matching) ;
        k = get_match(j, matching) ;
        VSP_ASSERT( IN_RANGE(j, 0, n) ) ;

        if(matchtype[i] != matchtype[j])
        {
            ERR("Matched vertices don't have the same matchtype") ;
        }

        switch (matchtype[i])
        {
            case matchtype_HEM :
                // Check that j is matched to i
                if( get_match(j, matching) != i )
                {
                    ERR("Error in HEM pair") ;
                }

                // Check the matchtype of j is HEM
                if(matchtype[j] != matchtype_HEM)
                {
                    ERR("Incorrect matchtype") ;
                }

                // Check that i, j are neighbors
                if( !Are_Neighbors(i, j, Ai, Ap) )
                {
                    ERR("Nodes in HEM pair are not neighbors") ;
                }
                break ;
            case matchtype_brotherly:
                // Check that j is matched to i
                if( get_match(j, matching) != i )
                {
                    ERR("Error in BM pair") ;
                }

                // Check the matchtype of j is BM
                if(matchtype[j] != matchtype_brotherly)
                {
                    ERR("Incorrect matchtype") ;
                }

                // Check that i, j are NOT neighbors
                if( Are_Neighbors(i, j, Ai, Ap) )
                {
                    ERR("Nodes in BM pair ARE neighbors") ;
                }
                break ;

            case matchtype_3way:
                // Check that i is in a 3-way match
                k = get_match(j, matching) ;
                VSP_ASSERT(k != i) ;
                VSP_ASSERT(k != j) ;
                if(get_match(k, matching) != i)
                {
                    ERR("Error in 3-way matching") ;
                }

                // Check that i neighbors at least one of j, k
                if(!Are_Neighbors(i, j, Ai, Ap) && !Are_Neighbors(i, k, Ai, Ap)) 
                {
                    ERR("Invalid 3-way match") ;
                }
                break ;
        }
    }

    Check_cn(n, cn, matching) ;
}


void Graph_Simple_Match
(
    Graph       *Graph,
    Matching        *M,
    VSPINT         lvl,
    Options   *options
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT *cn = &(M->cn) ;
    VSPINT *matching = M->match ;
    VSPINT *matchmap = M->matchmap ;
    VSPINT *invmatchmap = M->invmatchmap ;
    match_type *matchtype = M->matchtype ;

    VSPINT i, j, k, v, q, j_ew, k_ew, ck, cj, ci_coarse ;
    VSPLONG p, end ;

    // Graph info
    VSPINT   n   = Graph->n ;
    VSPLONG  nz  = Graph->nz ;
    VSPINT  *Ai  = Graph->Ai ;
    VSPINT  *Ae  = Graph->Ae ;
    VSPLONG *Ap  = Graph->Ap ;
    VSPINT  *c   = Graph->c ;
    VSPINT pow2l   = pow_ints(2, lvl) ;
    VSPINT pow2lp1 = pow_ints(2, lvl + 1) ;

    bool Ae_is_one = (Ae == NULL) ;

    VSPINT n_matches   = 0 ;
    VSPINT n_unmatched = n ;

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Heavy Edge and Weight Matching
    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    VSPINT *heap_sings      = malloc((n+1) * sizeof(*heap_sings)) ;
    VSPINT  n_heap_sings    = 0 ;

    VSPINT *lookup_heap     = malloc((n+1) * sizeof(*lookup_heap)) ;

    // The largest cost of any 2-way or 3-way match. Used to block singleton
    // matching if the cost of such a matching would exceed gamma_coarse.
    VSPINT gamma_coarse     = 0 ; 
    // bool exists_perfect_coarse_vtx = false ;

    for (i = 0; i < n; i++)
    {
        if(Is_Matched(i, matching)) continue ;
        if(Is_Singleton(i, Ap))
        {
            heap_sings[++n_heap_sings] = i ;
            lookup_heap[i] = n_heap_sings ;
            continue ;
        }

        // j holds the unmatched neighbor with heaviest edge and greatest
        // vertex cost
        j    = EMPTY ;
        j_ew = EMPTY ;
        cj   = EMPTY ;
        end  = Ap[i+1] ;

        for (p = Ap[i] ; p < end ; p++)
        {
            k = Ai[p];
            if(Is_Matched(k, matching)) continue ;

            k_ew = Ae_is_one ? 1 : Ae[p] ;
            if ( k_ew > j_ew || (k_ew == j_ew && c[k] > cj) )
            {
                j    = k ;
                j_ew = k_ew ;
                cj   = c[k] ;
            }
        }

        // Match i with j
        if (j != EMPTY)
        {
            Create_Match_2way
            (
                i, j, matching, invmatchmap, matchtype,
                matchmap, cn, matchtype_HEM
            ) ;
            n_matches++ ;
            n_unmatched -= 2 ;

            ci_coarse = c[i] + c[j] ;
            if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Match singletons
    /*------------------------------------------------------------------------*/
    VSP_ASSERT(n_heap_sings >= 0) ;

    // Can occur if the graph consists entirely of singletons
    if(gamma_coarse == 0) gamma_coarse = pow2lp1 ;
    Match_Singletons
    (
        heap_sings, lookup_heap, &n_heap_sings, c, pow2lp1, matching,
        matchmap, invmatchmap, matchtype, cn, &n_matches, &n_unmatched
    ) ;
    free(heap_sings) ;
    free(lookup_heap) ;
    heap_sings  = NULL ;
    lookup_heap = NULL ;

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Brotherly matching with pivots
    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    VSPINT node_pivot, ew_pivot, neigh, ew ;
    VSPINT perfect, light, lightest, c_lightest ;
    for(k = 0; k < n; k++)
    {
        if ( Is_Matched(k, matching) ) continue;

        node_pivot = EMPTY ;
        ew_pivot   = EMPTY ;

        // Find the neighbor of k with heaviest weight
        for (p = Ap[k]; p < Ap[k+1]; p++)
        {
            neigh = Ai[p] ;
            VSP_ASSERT(Is_Matched(neigh, matching)) ;

            ew = Ae_is_one ? 1 : Ae[p] ;
            if (ew > ew_pivot)
            {
                ew_pivot   = ew ;
                node_pivot = neigh;
            }
        }

        if(matchtype[node_pivot] == matchtype_3way)
        {
            printf("\nWARNING\n") ;
            fflush(stdout) ;
        }

        if(matchtype[node_pivot] != matchtype_HEM)
        {
            printf("\nWARNING2\n") ;
            fflush(stdout) ;
        }

        // Match neighbors of the pivot, prioritizing matches between perfect
        // vertices
        if (node_pivot != EMPTY)
        {
            perfect = EMPTY ;
            light   = EMPTY ;

            for (p = Ap[node_pivot] ; p < Ap[node_pivot + 1] ; p++)
            {
                neigh = Ai[p] ;
                if (Is_Matched(neigh, matching)) continue ;

                if(Is_Perfect_Vtx(neigh, c, pow2l))
                {
                    if(perfect != EMPTY)
                    {
                        Create_Match_2way
                        (
                            perfect, neigh, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        perfect = EMPTY ;
                        n_matches++ ;
                        n_unmatched -= 2 ;
                    }
                    else perfect = neigh ;
                }
                else
                {
                    if(light != EMPTY)
                    {
                        Create_Match_2way
                        (
                            light, neigh, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        n_matches++ ;
                        n_unmatched -= 2 ;
                        light = EMPTY ;
                    }
                    else light = neigh ;
                }
            }
        }
    }


    VSPINT n_coarse = n_matches + n_unmatched ;
    VSPFLOAT ratio_compression = ((VSPFLOAT) n_coarse) / ((VSPFLOAT) n) ;

    #ifdef DEBUG
        VSP_ASSERT(*cn == n_matches) ;
        VSPINT dbg_n_unmatched = 0 ;
        for(VSPINT ii = 0 ; ii < n ; ii++)
        {
            if(!Is_Matched(ii, matching)) dbg_n_unmatched++ ;
        }
        if(dbg_n_unmatched != n_unmatched)
        {
            ERR("Wrong count of unmatched nodes") ;
        }
    #endif


    if(ratio_compression > options->compression_threshold)
    {
        /*--------------------------------------------------------------------*/
        /*--------------------------------------------------------------------*/
        // Brotherly matching new
        /*--------------------------------------------------------------------*/
        /*--------------------------------------------------------------------*/

        VSPINT *selfmatch    = malloc(n * sizeof(*selfmatch)) ;
        VSPINT num_selfmatch = 0 ;

        // Loop through HEM pairs
        if((size_t) n > 2E9 || (size_t) n <= 0)
        {
            printf("\n\nInvalid n: %zu \n\n", (size_t) n) ;
        }
        bool *HEM_Pair_Visited = calloc(n, sizeof(bool)) ;
        if((size_t) n > 2E9 || (size_t) n <= 0)
        {
            printf("\n\nInvalid n: %zu \n\n", (size_t) n) ;
        }

        for(i  = 0 ; i < n ; i++)
        {
            if(matchtype[i] != matchtype_HEM || HEM_Pair_Visited[i] ) continue ;

            // {i, j} is a HEM pair 
            j = get_match(i, matching) ;

            // Mark {i, j} pair as visited
            HEM_Pair_Visited[i] = true ;
            HEM_Pair_Visited[j] = true ;

            // Initialize
            perfect    = EMPTY ;
            light      = EMPTY ;
            lightest   = EMPTY ;
            c_lightest = pow2l + 1 ;
            end = Ap[i+1] ;

            /*----------------------------------------------------------------*/
            // Visit Neighbors of i
            /*----------------------------------------------------------------*/
            for(p = Ap[i] ; p < end ; p++)
            {

                k  = Ai[p] ;
                if(Is_Matched(k, matching)) continue ;

                ck = c[k] ;
                if(ck < c_lightest)
                {
                    lightest = k ;
                    c_lightest = ck ;
                }

                if(ck == pow2l)
                {
                    if(perfect != EMPTY)
                    {
                        Create_Match_2way
                        (
                            perfect, k, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        ci_coarse = c[i] + c[j] ;
                        #ifdef DEBUG
                            if(c[perfect] + c[k] != pow2lp1)
                            {
                                ERR("Non-perfect vertex") ;
                            }
                        #endif
                        perfect = EMPTY ;
                    }
                    else perfect = k ;
                }
                else // ck != pow2l, i.e. k is not perfect
                {
                    VSP_ASSERT(ck < pow2l) ;
                    if(light != EMPTY)
                    {
                        Create_Match_2way
                        (
                            light, k, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        if(gamma_coarse < pow2l)
                        {
                            ci_coarse = c[i] + c[j] ;
                            if(ci_coarse > gamma_coarse) gamma_coarse=ci_coarse;
                        }
                        light = EMPTY ;
                    }
                    else light = k ;
                }
            }

            /*----------------------------------------------------------------*/
            // Visit Neighbors of j
            /*----------------------------------------------------------------*/
            end = Ap[j+1] ;
            for(p = Ap[j] ; p < end ; p++)
            {
                k  = Ai[p] ;
                if(Is_Matched(k, matching)) continue ;

                // This check is for the following situation: u is a neighbor of
                // both i and j. After visiting the neighbors of i, assume that 
                // perfect = u or light = u. Then at some point k will be u, and
                // if we don't stop this iteration, then k will be matched with
                // itself
                if(k == perfect || k == light) continue ;

                ck = c[k] ;
                if(ck < c_lightest)
                {
                    lightest = k ;
                    c_lightest = ck ;
                }

                if(ck == pow2l)
                {
                    if(perfect != EMPTY)
                    {
                        Create_Match_2way
                        (
                            perfect, k, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        perfect = EMPTY ;
                        ci_coarse = c[i] + c[j] ;
                    }
                    else perfect = k ;
                }
                else // ck == pow2l, i.e. k is not perfect
                {
                    VSP_ASSERT(ck < pow2l) ;
                    if(light != EMPTY)
                    {
                        Create_Match_2way
                        (
                            light, k, matching, invmatchmap, matchtype,
                            matchmap, cn, matchtype_brotherly
                        ) ;
                        light = EMPTY ;
                        if(gamma_coarse < pow2l)
                        {
                            ci_coarse = c[i] + c[j] ;
                            if(ci_coarse > gamma_coarse) gamma_coarse=ci_coarse;
                        }

                    }
                    else light = k ;
                }
            }

            /*----------------------------------------------------------------*/
            // Resolve unmatched neighbors of {i, j}
            /*----------------------------------------------------------------*/

            // If {i, j} has two unmatched neighbors remaining, match them
            if(perfect != EMPTY && light != EMPTY)
            {
                Create_Match_2way
                (
                    perfect, light, matching, invmatchmap, matchtype,
                    matchmap, cn, matchtype_brotherly
                ) ;

                ci_coarse = c[i] + c[j] ;
                if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;
                continue ;
            }

            // Here, at least one of perfect and light is non-EMPTY
            // If both are EMPTY, k is EMPTY. Otherwise, k is the non-empty one
            k = perfect + light + 1 ;

            if( k != EMPTY )
            {
                ck = c[k] ;

                // Ensure the remaining vertex has lightest cost
                if(ck > c_lightest)
                {
                    Steal_Match
                    (
                        lightest, k, matching, invmatchmap, matchtype, matchmap
                    ) ;

                    k  = lightest ;
                    ck = c_lightest ;
                }

                // If the 3-way match is light, do it
                if(c[i] + c[j] + c[k] <= pow2lp1)
                {
                    Create_Match_3way_HEM
                    (
                        i, j, k, matching, matchmap, invmatchmap, matchtype
                    ) ;

                    if(c[i] + c[j] + c[k] > gamma_coarse)
                    {
                        gamma_coarse = c[i] + c[j] + c[k] ;
                    }
                }
                else selfmatch[ num_selfmatch++ ] = k ;
            }
        }

        // Self-match any remaining vertices
        for(i = 0 ; i < num_selfmatch ; i++)
        {
            j = selfmatch[i] ;
            if( !Is_Matched(j, matching) )
            {
                Create_Match_2way
                (
                    j, j, matching, invmatchmap, matchtype,
                    matchmap, cn, matchtype_selfMatch
                ) ;
                if(c[j] > gamma_coarse) gamma_coarse = c[j] ;
            }
        }
        free(selfmatch) ;
        free(HEM_Pair_Visited) ;
        selfmatch        = NULL ;
        HEM_Pair_Visited = NULL ;
    }
    else
    {
        // Self-match any unmatched vertices
        for(i = 0 ; i < n ; i++)
        {
            if( !Is_Matched(i, matching) )
            {
                Create_Match_2way
                (
                    i, i, matching, invmatchmap, matchtype,
                    matchmap, cn, matchtype_selfMatch
                ) ;
            }
        }
    }

    /* ====================================================================== */
    // Cleanup
    /* ====================================================================== */
    #ifdef DEBUG
        Match_Check
        (
            Graph->n,
            M
        ) ;

        Additional_Match_Check
        (
            Graph,
            lvl,
            matching,
            matchmap,
            invmatchmap,
            matchtype,
            *cn
        ) ;
    #endif

    return ;
}


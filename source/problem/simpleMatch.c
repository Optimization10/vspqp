#include "shared_api.h"
#include "shared_functions.h"
#include "shared_macros.h"
#include "shared_types.h"
#include "coarsen_internal.h"
#include "ds_api.h"


static void Additional_Match_Check
(
    VSP_Graph             *Graph,
    VSPINT                   lvl,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    VSP_match_type    *matchtype
) ;

static void Create_Match_2way
(
    VSPINT                       i,
    VSPINT                       j,
    VSPINT*               matching,
    VSPINT*            invmatchmap,
    VSP_match_type*  matchtype_arr,
    VSPINT*               matchmap,
    VSPINT*                 cn_ptr,
    VSP_match_type           mtype
)
{
    SET_FUNCTION_COUNTER;
    #ifdef VSP_DEBUG
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
    VSP_match_type*    matchtype_arr
)
{

    SET_FUNCTION_COUNTER;
    #ifdef VSP_DEBUG
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
    VSP_match_type*    matchtype_arr,
    VSPINT*            cn_ptr
)
{

    SET_FUNCTION_COUNTER;
    #ifdef VSP_DEBUG
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
    VSP_match_type*      matchtype,
    VSPINT*               matchmap
)
{
    VSPINT i = VSP_get_match(old, matching) ;
 
    // Update matching
    matching[old] = 0 ;
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

extern struct timespec  G_START_TIME, G_END_TIME ;
extern long double G_HEM_TIME, G_BM_TIME ;


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
    VSP_match_type    *matchtype,
    VSPINT                   *cn
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
        }
        else if(j != EMPTY && c[i] + c[j] <= gamma_coarse)
        {
            // Match i with j
            Create_Match_2way
            (
                i, j, matching, invmatchmap, matchtype,
                matchmap, cn, matchtype_singleton
            ) ;
        }
        else
        {
            // Match i with itself
            Create_Match_2way
            (
                i, i, matching, invmatchmap, matchtype,
                matchmap, cn, matchtype_singleton
            ) ;

            if(j != EMPTY)
            {
                // Match j with itself
                Create_Match_2way
                (
                    j, j, matching, invmatchmap, matchtype,
                    matchmap, cn, matchtype_singleton
                ) ;
            }
        }
    }
}



void Graph_Simple_Match
(
    VSP_Graph              *Graph,
    VSPINT                    lvl,
    VSPINT              *matching,
    VSPINT              *matchmap,
    VSPINT           *invmatchmap,
    VSP_match_type     *matchtype,
    VSPINT                    *cn
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, j, k, q, j_ew, k_ew, ck, cj, ci_coarse ;
    VSPLONG p, end ;

    // Graph info
    VSPINT   n   = Graph->n ;
    VSPLONG  nz  = Graph->nz ;
    VSPINT  *Ai  = Graph->Ai ;
    VSPINT  *Ae  = Graph->Ae ;
    VSPLONG *Ap  = Graph->Ap ;
    VSPINT  *c   = Graph->c ;
    VSPINT pow2l = pow_ints(2, lvl + 1) ;


    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Heavy Edge and Singleton Matching
    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    VSPINT *heap_sings      = malloc((n+1) * sizeof(*heap_sings)) ;
    VSPINT  n_heap_sings    = 0 ;

    VSPINT *lookup_heap      = malloc((n+1) * sizeof(*heap_sings)) ;

    // The largest cost of any 2-way or 3-way match. Used to block singleton
    // matching if the cost of such a matching would exceed gamma_coarse.
    VSPINT gamma_coarse     = 0 ; 
    // bool exists_perfect_coarse_vtx = false ;

    clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;
    for (i = 0; i < n; i++)
    {
        if(Is_Matched(i, matching)) continue ;
        if(Is_Singleton(i, Ap))
        {
            heap_sings[++n_heap_sings] = i ;
            lookup_heap[i] = n_heap_sings ;
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

            k_ew = Ae[p] ;
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
            ci_coarse = c[i] + c[j] ;
            if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;
        }
    }

    clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
    G_HEM_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;

    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    // Brotherly, Light 3-way, and Self Matching
    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;

    // Unmatched perfect and non-perfect neighbors
    VSPINT perfect, light, lightest, c_lightest ;

    // Hold vertices which may need to be self-matched
    VSPINT *selfmatch    = malloc(n * sizeof(*selfmatch)) ;
    VSPINT num_selfmatch = 0 ;

    // Loop through HEM pairs
    bool *HEM_Pair_Visited = calloc(n, sizeof(*HEM_Pair_Visited)) ;
    for(i  = 0 ; i < n ; i++)
    {
        if(matchtype[i] != matchtype_HEM || HEM_Pair_Visited[i] ) continue ;

        // {i, j} is a HEM pair 
        j = VSP_get_match(i, matching) ;

        // Mark {i, j} pair as visited
        HEM_Pair_Visited[i] = true ;
        HEM_Pair_Visited[j] = true ;

        // Initialize
        perfect    = EMPTY ;
        light      = EMPTY ;
        lightest   = EMPTY ;
        c_lightest = pow2l + 1 ;
        end = Ap[i+1] ;

        /*--------------------------------------------------------------------*/
        // Visit Neighbors of i
        /*--------------------------------------------------------------------*/
        for(p = Ap[i] ; p < end ; p++)
        {

            k  = Ai[p] ;
            if(Is_Matched(k, matching)) continue ;

            ck = c[k] ;
            if(ck < c_lightest) lightest = k ;

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
                    #ifdef VSP_DEBUG
                        if(c[perfect] + c[k] != pow2l)
                        {
                            ERR("Non-perfect vertex") ;
                        }
                    #endif
                    perfect = EMPTY ;
                }
                else perfect = k ;
            }
            else // ck == pow2l, i.e. k is not perfect
            {
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
                        if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;
                    }
                    light = EMPTY ;
                }
                else light = k ;
            }
        }

        /*--------------------------------------------------------------------*/
        // Visit Neighbors of j
        /*--------------------------------------------------------------------*/
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
            // I wonder how expe
            if(k == perfect || k == light) continue ;

            ck = c[k] ;
            if(ck < c_lightest) lightest = k ;

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
                        if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;
                    }

                }
                else light = k ;
            }
        }

        /*--------------------------------------------------------------------*/
        // Resolve unmatched neighbors of {i, j}
        /*--------------------------------------------------------------------*/

        // If {i, j} has two unmatched neighbors remaining, match them together
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

                // Check if gamma_coarse increased due to stolen match
                ci_coarse = c[k] + c[ VSP_get_match(k, matching) ] ;
                if(ci_coarse > gamma_coarse) gamma_coarse = ci_coarse ;

                k  = lightest ;
                ck = c_lightest ;
            }

            // If the 3-way match is light, do it
            if(c[i] + c[j] + c[k] <= pow2l)
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

    /* ====================================================================== */
    // Match singletons
    /* ====================================================================== */
    VSP_ASSERT(n_heap_sings >= 0) ;
    Match_Singletons
    (
        heap_sings, lookup_heap, &n_heap_sings, c, gamma_coarse, matching,
        matchmap, invmatchmap, matchtype, cn
    ) ;


    /* ====================================================================== */
    // Cleanup
    /* ====================================================================== */
    free(HEM_Pair_Visited) ;
    free(selfmatch) ;
    free(heap_sings) ;
    free(lookup_heap) ;

    clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
    G_BM_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;

    #ifdef VSP_DEBUG
        Match_Check
        (
            n,
            *cn,
            matching,
            matchmap,
            invmatchmap,
            matchtype
        ) ;

        Additional_Match_Check
        (
            Graph,
            lvl,
            matching,
            matchmap,
            invmatchmap,
            matchtype
        ) ;
    #endif
}


// Check the following:
// A vertex and its match(es) all have the same matchtype

// If i has matchtype HEM, then it's matched to precisely one other vertex,
// and this vertex is a neihbor

// If i has matchtype brotherly, then it is matched to precisely one other
// vertex and this vertex is NOT a neighbor

// If i is in a 3-way match, then it is matched to precisely 2 other vertices,
// and it neighbors at least one of them

static void Additional_Match_Check
(
    VSP_Graph             *Graph,
    VSPINT                   lvl,
    VSPINT             *matching,
    VSPINT             *matchmap,
    VSPINT          *invmatchmap,
    VSP_match_type    *matchtype
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
        j = VSP_get_match(i, matching) ;
        k = VSP_get_match(j, matching) ;
        VSP_ASSERT( IN_RANGE(j, 0, n) ) ;

        if(matchtype[i] != matchtype[j])
        {
            ERR("Matched vertices don't have the same matchtype") ;
        }

        switch (matchtype[i])
        {
            case matchtype_HEM :
                // Check that j is matched to i
                if( VSP_get_match(j, matching) != i )
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
                if( VSP_get_match(j, matching) != i )
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
                k = VSP_get_match(j, matching) ;
                VSP_ASSERT(k != i) ;
                VSP_ASSERT(k != j) ;
                if(VSP_get_match(k, matching) != i)
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
}
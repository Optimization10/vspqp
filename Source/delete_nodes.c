#include "delete_nodes.h"
#include "shared_api.h"
#include "shared_functions.h"
#include "shared_types.h"

// Return NULL if the input graph is all singletons
VSP_Graph *Compress_Graph
(
    VSP_Graph *GraphOrg,
    VSPINT        *keep,   // keep[i] = 1 if to keep i, 0 if to not keep i
    VSPINT    *newIndex    // Output
)
{
    VSPINT j, i, compressedVtx,
    nOrg = GraphOrg->n ;
    VSPLONG k, l, m,
    nzOrg = GraphOrg->nz ;

    // Construct reduced graph to be returned
    VSP_Graph *ret = Construct_Graph(nOrg, nzOrg) ;    
    
    // Arrays containing original graph
    VSPLONG *ApOrg = GraphOrg->Ap ;
    VSPINT  *AiOrg = GraphOrg->Ai ;
    VSPINT  *AeOrg = GraphOrg->Ae ;
    VSPINT   *cOrg = GraphOrg->c ;

    // Arrays to contain reduced graph
    VSPLONG *ApRed = ret->Ap ;
    VSPINT  *AiRed = ret->Ai ;
    VSPINT  *AeRed = ret->Ae ;
    VSPINT  *cRed  = ret->c ;

    VSPINT minCostRed = GraphOrg->total_vtx_cost ; // Upper bound for minCostRed
    VSPINT maxCostRed = 0 ;
    VSPINT totalVtxCostRed = 0 ;

    // Get the reduced cost vector, and reduced min, max, and total vertex cost
    // Initialize newIndex to all EMPTY's
    j = 0 ;
    for (i = 0; i < nOrg; i++)
    {
        newIndex [i] = EMPTY ;
        if ( keep [i] )
        {
            newIndex [i] = j ; /* j = compressed index */
            cRed [j] = cOrg [i] ;
            if(cRed[j] > maxCostRed) maxCostRed = cRed[j] ;
            if(cRed[j] < minCostRed) minCostRed = cRed[j] ;
            totalVtxCostRed += cRed[j] ;
            j++ ;
        }
    }
    ret->min_vtx_cost   = minCostRed ;
    ret->max_vtx_cost   = maxCostRed ;
    ret->total_vtx_cost = totalVtxCostRed ;

    // Build the reduced Ai, Ap, Ae
    l = 0 ;
    compressedVtx = 0 ;
    ApRed[0] = 0 ;
    for (j = 1 ; j <= nOrg ; j++)
    {
        if ( ! keep[j-1] ) continue ;

        m = ApOrg[j] ; /* end of column j-1 */
        for (k = ApOrg[j-1]; k < m; k++)
        {
            i = AiOrg[k] ;
            if (newIndex[i] != EMPTY)
            {
                AiRed[l] = newIndex[i] ;
                AeRed[l] = AeOrg[k] ;
                l++ ;
            }
        }
        compressedVtx++ ;
        ApRed[compressedVtx] = l ;
    }

    // Set the n and nz for the compressed graph
    ret->n  = compressedVtx ;
    ret->nz = l ;

    if(ret->n == 0)
    {
        Destroy_Graph(ret) ;
        ret = NULL ;
    }

    return ret ;
}


VSP_state *Compress_Binary_State
(
    VSP_Graph        *compGraph,
    VSPINT            *newIndex,
    VSP_state             *Zorg
)
{
    if(Zorg->point->ind_frac != EMPTY)
    {
        printf("\nERROR: Must compress binary point") ;
        return NULL ;
    }

    VSPINT norg = Zorg->point->n ;
    VSPINT ncmp = compGraph->n ;

    VSP_state *Zcmpr = VSP_construct_state(ncmp) ;

    VSPINT i, j,
    *zorg     = Zorg->point->int_part,
    *zcmpr    = Zcmpr->point->int_part,
    *ccmpr    = compGraph->c,
    zcostcmpr = 0 ;

    /*------------------------------------------------------------------------*/
    // Compress z and cost of z
    /*------------------------------------------------------------------------*/
    for(i = 0 ; i < norg ; i++)
    {
        j = newIndex[i] ;
        if(j != EMPTY)
        {
            zcmpr[j] = zorg[i] ;
            zcostcmpr += zcmpr[j] * ccmpr[j] ; 
        }
    }
    Zcmpr->weight = zcostcmpr ;

    /*------------------------------------------------------------------------*/
    // Compress Az
    /*------------------------------------------------------------------------*/
    VSPINT  *Azcmp = Zcmpr->Az->int_part  ;
    VSPINT  *Aicmp = compGraph->Ai ;
    VSPLONG *Apcmp = compGraph->Ap ;

    VSP_mat_vec_mult(zcmpr, Azcmp, Aicmp, Apcmp, ncmp) ;

    return Zcmpr ;
}




// Given a state for the compressed problem, return the corresponding state
// for the uncompressed problem where the vertices which were not kept are in
// the separator 
VSP_state *Uncompress_State
(
    VSP_Graph        *compGraph,
    VSPINT            *newIndex,
    VSP_state             *Zcmp,
    VSPINT               *AiOrg,
    VSPLONG *ApOrg,
    VSPINT                 nOrg
)
{
    VSPINT i, j, *zOrg,
    *zCmp = Zcmp->point->int_part ;

    /*------------------------------------------------------------------------*/
    // Allocate memory for uncompressed state
    /*------------------------------------------------------------------------*/
    VSP_state *ret = VSP_construct_state(nOrg) ;
    zOrg = ret->point->int_part ;
    

    /*------------------------------------------------------------------------*/
    // Uncompress z and weight
    /*------------------------------------------------------------------------*/
    // Set the integer part. No need to set uncompressed entries to 0 explicitly
    // since the array is calloc'ed

    VSPINT cmpIndFrac = Zcmp->point->ind_frac ;
    for(i = 0 ; i < nOrg ; i++)
    {
        if( (j = newIndex[i]) != EMPTY)
        {
            // i in the uncompressed problem <---> j in the compressed problem
            zOrg[i] = zCmp[j] ;

            // Set the fractional part
            if(cmpIndFrac != EMPTY)
            {
                if(j == cmpIndFrac)
                {
                    ret->point->ind_frac = i ;
                    ret->point->val_frac = Zcmp->point->val_frac ;
                }
            }
        }
    }

    // Uncompress the weight
    ret->weight = Zcmp->weight ;

    /*------------------------------------------------------------------------*/
    // Uncompress Az
    /*------------------------------------------------------------------------*/
    Compute_Az(ret->point, ret->Az, AiOrg, ApOrg, nOrg) ;

    return ret ;
}
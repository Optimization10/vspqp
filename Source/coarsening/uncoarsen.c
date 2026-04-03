#include "core.h"
#include "shared_api.h"
#include "coarsen_internal.h" // IWYU pragma: keep
#include "shared_functions.h"
#include "shared_macros.h"
#include <stdio.h>



static void CheckForFracEntryAndRefineIt
(
    VSP_state        *Zcrs,
    VSP_state       *Zfine,
    Problem  *probFine
) ;


static void RefineOneVar
(
    Problem     *fineGraph,
    VSP_state           *fineZ,
    Problem   *coarseGraph, 
    VSP_state         *coarseZ
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;

        /* Check input */
        VSP_check_state(coarseZ, coarseGraph) ;
    #endif

    VSPINT i, fn, fu, fv, cv, fvWgt, *fc, *cz, *fz, *matchMap ;
    VSPFLOAT fWgtz ;

    fn    = fineGraph ->Graph->n ;
    fc    = fineGraph ->Graph->c ;
    cz    = coarseZ   ->point->int_part ;
    fz    = fineZ     ->point->int_part ;
    fu    = fineGraph->u ;
    fWgtz = 0 ;
    matchMap = fineGraph->matchmap ;

    /*------------------------------------------------------------------------*/
    // Set all entries of fine Az to 0
    /*------------------------------------------------------------------------*/
    VSPINT *fAz = fineZ->Az->int_part ;
    for(i = 0 ; i < fn ; i++) fAz[i] = 0 ;
    fineZ->Az->ind_frac_z = EMPTY ;
    fineZ->Az->val_frac_z = 0 ;


    /*------------------------------------------------------------------------*/
    /*   Refine the Integer Part */
    /*------------------------------------------------------------------------*/
    for( fv = 0 ; fv < fn ; fv++ )
    {
        cv = matchMap[fv];
        if( cz[ cv ] ) /* The VSPINT part of the coarse vertex is 1 */
        {
            fvWgt = fc[ fv ] ;
            if( fWgtz + fvWgt < fu ) /* fv fits into z with extra room */
            {
                fz[ fv ] = 1 ;
                fWgtz += fvWgt ;
                VSP_Az_pm_col_of_A(fAz, fv, PLUS, fineGraph) ;
                continue ;
            }
            else if ( fWgtz + fvWgt == fu  ) /* fv fits, but no extra room */
            {
                fz[ fv ] = 1 ;
                fWgtz = fu ;
                VSP_Az_pm_col_of_A(fAz, fv, PLUS, fineGraph) ;
                break ;
            }
            else /* fv partially fits into z */
            {
                fWgtz    = fu ;
                break ;
            }
        }
        else /* The VSPINT part of the coarse vertex is 0 */
        {
            fz[fv] = 0 ;
        }
    }

    /* Set any remaining fine vertices to 0 */
    for( fv++ ; fv < fn ; fv++ )
    {
        fz[fv] = 0 ;
    }

    fineZ->weight = fWgtz ;

    /*------------------------------------------------------------------------*/
    /*   Refine the Frac Entry and Update Az */
    /*------------------------------------------------------------------------*/
    CheckForFracEntryAndRefineIt(coarseZ, fineZ, fineGraph) ;

    #ifdef VSP_DEBUG
        VSP_check_state(fineZ, fineGraph) ;
    #endif

    return ;
}


/*  This function assumes that
    1) Zfine contains the refined integer part. In particular, the entries in
       Zfine corresponding to the fractional entry in Zcrs (if it even exists)
       are 0

    This function does NOT touch Az.

    This function goes through the refined vertices and sets to 1 as many
    of them as possible. The function does NOT take into account the weight
    of the fractional entry in the coarse point. If the combined weight of
    the fine vertixes exceeds the weight of the fractional component, so be it
    (this can happen if the u for the fine problem is larger than the u for the
    coarse problem)

    */
static void CheckForFracEntryAndRefineIt
(
    
    VSP_state        *Zcrs,
    VSP_state       *Zfine,
    Problem  *probFine
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT ufine = probFine->u ;
    VSPINT crsIndFrac = Zcrs->point->ind_frac ;
    VSPINT *cfine = probFine->Graph->c ;
    VSPINT *zfineInt  = Zfine->point->int_part ;

    /* If there is no fractional entry in the coarse point, or if there is no
       room in the fine point, there is nothing to do */
    if(crsIndFrac == EMPTY || Zfine->weight == ufine)
    {
        return ;
    }
    VSP_ASSERT( Zfine->point->ind_frac == EMPTY ) ;

    VSPINT firstFineVert = probFine->invmatchmap[crsIndFrac] ;
    VSPINT fineVert = firstFineVert ;
    while( 1 )
    {
        /* Set z[fineVert] */
        if( Zfine->weight + cfine[fineVert] < ufine )
        {
            zfineInt[fineVert] = 1 ;
            Zfine->weight += cfine[fineVert] ;
        }
        else if( Zfine->weight + cfine[fineVert] == ufine )
        {
            zfineInt[fineVert] = 1 ;
            Zfine->weight = ufine ;

            Zfine->point->ind_frac = EMPTY ;
            Zfine->Az->ind_frac_z  = EMPTY ;
            Zfine->point->val_frac = 0 ;
            Zfine->Az->val_frac_z  = 0 ;
            break ;   
        }
        else
        {
            VSPFLOAT valFrac = (ufine - Zfine->weight) / ( (VSPFLOAT) cfine[fineVert] ) ;
            Zfine->point->ind_frac = fineVert ;
            Zfine->point->val_frac = valFrac ;
            Zfine->weight = ufine ;
            break ;
        }

        fineVert = VSP_get_match(fineVert, probFine->matching) ;
        if(fineVert == firstFineVert)
        {
            break ;
        }
    }
    return ;
}



void Uncoarsen_Point
(
    Problem    *fine_graph,
    VSP_state          *fine_X,
    VSP_state          *fine_Y,
    VSP_state        *coarse_X,
    VSP_state        *coarse_Y
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    Problem *coarse_graph = fine_graph->child_problem ;
    RefineOneVar(fine_graph, fine_X, coarse_graph, coarse_X) ;
    RefineOneVar(fine_graph, fine_Y, coarse_graph, coarse_Y) ;

    #ifdef VSP_DEBUG
        VSP_check_states(fine_X, fine_Y, fine_graph) ;

        if(fine_graph->u >= coarse_graph->u)
        {
            VSP_ASSERT(fine_X->weight >= coarse_X->weight) ;
            VSP_ASSERT(fine_Y->weight >= coarse_Y->weight) ;
        }
    #endif 

    VSP_destruct_state(coarse_X) ;
    VSP_destruct_state(coarse_Y) ;
    return ;
}
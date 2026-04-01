#include "core.h"
/*----------------------------------------------------------------------------*/
// Basic Flips Strategies module

// Purpose: The purpose of this module is to implement the component strategies
// of FM Flips v8, such as doing all beneficial flips
/*----------------------------------------------------------------------------*/

#include "basic_flip_strats.h"
#include "flips_context.h"
#include "ds_api.h"

// Do all beneficial flips in either var 0, var 1, or both (depending on
// whichVar). If withLocking is true, then lock the flips we do here.
void BasicFlipStrats_FeasBeneficialFlips
(
    FlipsCtx                  *FC,
    VSPINT              whichVars,
    bool              withLocking,
    VSPFLOAT                gamma
)
{
    IFDBG(FlipsCtx_Check(FC)) ;

    VSPINT nflips = 0 ;
    while( true )
    {
        Flip flip = FlipsCtx_Select_Flip(FC, whichVars, FEAS_BEN_FLIP, gamma) ;

        if(FlipsCtx_Is_Nonempty_Flip(&flip))
        {
            FlipsCtx_Apply_Flip(FC, &flip, withLocking) ;
            nflips++ ;
        }
        else break ;
    }
    IFDBG(FlipsCtx_Check(FC)) ;
}

// Do upflips until the target weight is just about to be surpassed
// ASSUMED:
// - whichVar = 0 or 1, but not 2
// - All downflips (so all zi = 1) are locked in whichVar
//
// NOT ASSUMED:
// - Not assumed that all upflips are unlocked
// void BasicFlipStrats_FeasUpflipsUntilTargetWgt
// (
//     FlipsCtx                  *FC,
//     VSPINT               whichVar, // Can NOT be both variables
//     bool              withLocking,
//     Flip_Selector_Mode       mode,
//     VSPINT              targetWgt
// )
// {
//     SET_FUNCTION_COUNTER ;
//     Flip flip ;
//     VSPINT *wgts = FlipsCtx_Get_Wgts(FC) ;

//     /*------------------------------------------------------------------------*/
//     // Check inputs
//     /*------------------------------------------------------------------------*/
//     VSP_ASSERT(whichVar == 0 || whichVar == 1) ;

//     #ifdef VSP_DEBUG
//         // // Check that all downflips are locked
//         // for(i = 0 ; i < n ; i++)
//         // {
//         //     if(vars[whichVar][i] && !locked[whichVar][i])
//         //     {
//         //         ERR("Downflip not locked") ;
//         //         break ;
//         //     }
//         // }
//     #endif

//     // At the start of the loop, this is the number of flips done so far
//     VSPINT nflips = 0 ;
//     while( true )
//     {
//         // Select flip
//         flip = FlipsCtx_Select_Flip(FC, whichVar , mode) ;

//         if( FlipsCtx_Is_Feasible_Flip(FC, &flip, targetWgt) )
//         {
//             FlipsCtx_Apply_Flip(FC, &flip, withLocking) ;
//             nflips++ ;
//         }
//         else break ;
//     }
//     return ;
// }
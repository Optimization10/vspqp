#include "core.h"
#include "greedy_flips_internal.h"
#include "data_structures.h"
#include "utils.h"

typedef enum
{
    IMPROVE_OBJ_ONLY_v5,
    IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5,
    GAIN_GEQ_0_v5
} Expected_Flip_Effect_v5 ;



// Flip Selector Check
void FlipSelector_Check_v5
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    char                  givenVarFlip,
    VSPINT                givenIndFlip,
    Problem               *problem,
    Expected_Flip_Effect_v5 expFlipEffect
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    // Validate input
    if(givenIndFlip != EMPTY)
    {
        if(givenVarFlip != 'x' && givenVarFlip != 'y')
        {
            ERR("Mismatch between given ind flip and given var flip") ;
        }
    }
    else
    {
        if(givenVarFlip != '\0')
        {
            ERR("Mismatch between given ind flip and given var flip") ;
        }
    }

    VSPINT i, k, xi, yi, Axi, Ayi, ci, bestGain = 0 ;
    VSPINT n = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c = problem->graph->c ;
    VSPINT  u = problem->u ;

    /*------------------------------------------------------------------------*/
    // Load the given flip into a Flip data structure of easier debugging
    /*------------------------------------------------------------------------*/
    VSPINT givenGain = 0 ;
    if( givenVarFlip == 'x' )
    {
        givenGain = GAIN_FLIP( DIR_FLIP(x[givenIndFlip]), c[givenIndFlip], gamma, Ay[givenIndFlip] ) ;
    }
    else if( givenVarFlip == 'y' )
    {
        givenGain = GAIN_FLIP( DIR_FLIP(y[givenIndFlip]), c[givenIndFlip], gamma, Ax[givenIndFlip] ) ;
    }

    Flip givenFlip = {0} ;
    givenFlip.gamma = gamma ;
    givenFlip.i    = givenIndFlip ;
    givenFlip.var  = givenVarFlip ;
    givenFlip.gain = givenGain ;
    givenFlip.u = u ;
    if(givenVarFlip == 'x')
    {
        VSP_ASSERT(givenIndFlip != EMPTY) ;
        givenFlip.wgt = wgtx ;
        givenFlip.zi   = x[givenIndFlip] ;
        givenFlip.Awi  = Ay[givenIndFlip] ;
        givenFlip.dir  = 1 - 2 * x[givenIndFlip] ;
        givenFlip.ci   = c[givenIndFlip] ;
        if(x[givenIndFlip] == 0)
        {
            givenFlip.dir = 1 ;
        }
        else
        {
            givenFlip.dir = -1 ;
        }
    }
    else if(givenVarFlip == 'y')
    {
        VSP_ASSERT(givenIndFlip != EMPTY) ;
        givenFlip.wgt = wgty ;
        givenFlip.zi   = y[givenIndFlip] ;
        givenFlip.Awi  = Ax[givenIndFlip] ;
        givenFlip.dir  = 1 - 2 * y[givenIndFlip] ;
        givenFlip.ci   = c[givenIndFlip] ;
    }
    else
    {
        givenVarFlip = '\0' ;
        givenFlip.i  = EMPTY ;
    }


    // Check for feasibility
    if(givenFlip.i != EMPTY)
    {
        if(givenFlip.dir == 1 && !FEASIBLE_UPFLIP(givenFlip.wgt, givenFlip.ci, givenFlip.u) )
        {
            ERR("Infeasible flip") ;
        }
    } 


    /*------------------------------------------------------------------------*/
    // Determine the best feasible flip. It is recorded and stored in trueFlip
    // regardless of the gain
    // existsAFeasNeutralDnflip - if the best gain is 0, trueFlip may be recorded
    //                            as a neutral upflip. But we may want to know if
    //                            there exists a neutral downflip. This var is true
    //                            or false regardless if the gain is positive or w/e
    /*------------------------------------------------------------------------*/


    Flip trueFlip = {0} ;
    trueFlip.var  = '\0' ;
    trueFlip.i    = EMPTY ;
    trueFlip.gamma = gamma ;
    trueFlip.u = u ;

    VSPINT gain, maxFeasGain = - n * gamma ;

    VSPINT feasNeutralDnFlipWitness = EMPTY ;
    for(i = 0 ; i < n ; i++)
    {
        xi  = x[i] ;
        yi  = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;
        ci  = c[i] ;
        if(xi)
        {
            gain = GAIN_FLIP(DOWN, ci, gamma, Ayi) ;
            if(gain > maxFeasGain)
            {
                maxFeasGain = gain ;
                trueFlip.Awi = Ayi ;
                trueFlip.ci = ci ;
                trueFlip.dir = DOWN ;
                trueFlip.gain = gain ;
                trueFlip.i = i ;
                trueFlip.var = 'x' ;
                trueFlip.wgt = wgtx ;
                trueFlip.zi = 1 ;
            }
            if(gain == 0 && wgtx > ci) feasNeutralDnFlipWitness = i ;
        }
        else // xi = 0
        {
            gain = GAIN_FLIP(UP, ci, gamma, Ayi) ;
            if(gain > maxFeasGain && FEASIBLE_UPFLIP(wgtx, ci, u))
            {
                maxFeasGain = gain ;
                trueFlip.Awi = Ayi ;
                trueFlip.ci = ci ;
                trueFlip.dir = UP ;
                trueFlip.gain = gain ;
                trueFlip.i = i ;
                trueFlip.var = 'x' ;
                trueFlip.wgt = wgtx ;
                trueFlip.zi = 0 ;
            }
        }

        if(yi)
        {
            gain = GAIN_FLIP(DOWN, ci, gamma, Axi) ;
            if(gain > maxFeasGain)
            {
                maxFeasGain = gain ;
                trueFlip.Awi = Axi ;
                trueFlip.ci = ci ;
                trueFlip.dir = DOWN ;
                trueFlip.gain = gain ;
                trueFlip.i = i ;
                trueFlip.var = 'y' ;
                trueFlip.wgt = wgty ;
                trueFlip.zi = 1 ;
            }
            if(gain == 0 && wgty > ci) feasNeutralDnFlipWitness = i ;
        }
        else // yi = 0
        {
            gain = GAIN_FLIP(UP, ci, gamma, Axi) ;
            if(gain > maxFeasGain && FEASIBLE_UPFLIP(wgty, ci, u))
            {
                maxFeasGain = gain ;
                trueFlip.Awi = Axi ;
                trueFlip.ci = ci ;
                trueFlip.dir = UP ;
                trueFlip.gain = gain ;
                trueFlip.i = i ;
                trueFlip.var = 'y' ;
                trueFlip.wgt = wgty ;
                trueFlip.zi = 0 ;
            }
        }
    }

    if(trueFlip.i == EMPTY)
    {
        printf("\n\n\nWARNING: VERY WEIRD CASE, NO FEASIBLE FLIP\n\n\n") ;
        ERR("\n\n\nWARNING: VERY WEIRD CASE, NO FEASIBLE FLIP\n\n\n") ;
        if(givenFlip.i != EMPTY)
        {
            ERR("") ;
        }
    }

    switch (expFlipEffect)
    {
    case IMPROVE_OBJ_ONLY_v5 :
        /*--------------------------------------------------------------------*/
        // If the true flip has a positive gain, check that the given flip has
        // the same gain. Otherwise, check that the given flip is empty
        /*--------------------------------------------------------------------*/
        if(trueFlip.gain > 0)
        {
            if((givenFlip.gain != trueFlip.gain) || givenFlip.i == EMPTY)
            {
                ERR("Gains mismatch") ;
            }
        }
        else
        {
            if(givenFlip.i != EMPTY)
            {
                ERR("Feasible flip that improves the obj doesn't actually exist") ;
            }
        }


        // Sanity check
        if(givenFlip.i != EMPTY && givenFlip.gain <= 0)
        {
            ERR("Bad flip selected, and code above is wrong") ;
        }

        break ;
    case IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5 :
        /*--------------------------------------------------------------------*/
        // Check the following:
        // If the true flip has positive gain, the given flip should too
        // If the true flip has gain 0 and there exists a feasible neutral
        // downflip, then given flip should be a neutral downflip
        // Otherwise, the given flip should be empty
        /*--------------------------------------------------------------------*/    
        if(trueFlip.gain > 0)
        {
            if( (trueFlip.gain != givenFlip.gain) || givenFlip.i == EMPTY)
            {
                ERR("Gains mismatch") ;
            }
        }
        else if(trueFlip.gain == 0 && feasNeutralDnFlipWitness != EMPTY)
        {
            if(givenFlip.i == EMPTY || givenFlip.gain != 0 || givenFlip.dir != -1)
            {
                ERR("Problem with given flip") ;
            }
        }
        else
        {
            if(givenFlip.i != EMPTY)
            {
                ERR("Problem with given flip") ;
            }
        }

        // Sanity check
        if(givenFlip.gain < 0)
        {
            ERR("Damaging flip selected, and code above is wrong") ;
        }
        break ;
    
    case GAIN_GEQ_0_v5 :
    {
       printf("\n\nUnsupported case\n") ;
       ERR("Unsupported case") ;
    }

    break ;
    default:
        ERR("Unknown expFlipEffect value") ;
        break;
    }
}


// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
//
// Fix vs. original v1: the gain-recording block has been split. Originally,
// `if(gain > maxFeasPosGain)` was doing both the "beneficial threshold"
// filter and the "new max" check at once, which caused any flip whose gain
// tied the running max to be silently dropped from the pool. Now the two
// jobs are separated:
//   - `if(gain > threshold)` records the flip (always, when beneficial)
//   - `if(gain > maxFeasPosGain)` updates the max (only on strict improvement)
//
// Pool building also now guards on `numFeasPosGain[flip] == 0` so that
// unrecorded flips can't accidentally tie when maxFeasPosGain == 0 in
// non-strict mode.
void SelectRandomOptimalBeneficialFlip
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,
    bool deep_search,
    bool strict,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list, // only used in debugging

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list, // only used in debugging

    Problem       *problem,
    char           *retFlipVar, // '\0' if no flip
    VSPINT         *retIndFlip  // EMPTY if no flip
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        Check_All_Data
        (
            x, Ax, wgtx, y, Ay, wgty,
            x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
            y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
            problem
        ) ;
    #endif

    VSPINT u      = problem->u ;
    VSPINT gamma  = problem->graph->max_vtx_cost ;
    // Get a random number
    // VSPINT ranNum = rand() ;

    /*------------------------------------------------------------------------*/
    // Collect the following data
    // feasPosGains[k]   - gain of flip k (valid iff numFeasPosGain[k] > 0)
    // numFeasPosGain[k] - number of indices achieving that gain (0 means
    //                     flip k was not recorded, doubles as a "present" flag)
    // maxFeasPosGain    - max over recorded gains
    //
    // A flip is beneficial iff
    //   strict     : gain > 0
    //   non-strict : gain > -1   (gain >= 0)
    /*------------------------------------------------------------------------*/
    VSPINT cval, gain, numIdx ;
    VSPINT threshold      = strict ? 0 : -1 ;
    VSPINT maxFeasPosGain = threshold ;
    VSPINT feasPosGains  [4] = {0} ;
    VSPINT numFeasPosGain[4] = {0} ;

    // x upflip
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, deep_search, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > threshold)
        {
            feasPosGains[XUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(x0Ay0_groups, cval, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }

    // y upflip
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, deep_search, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > threshold)
        {
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(y0Ax0_groups, cval, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }

    // x downflip
    IdxGroupsByC_See_Best_Feasible_Cval(x1Ay2_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > threshold)
        {
            feasPosGains[XDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x1Ay2_groups, &numIdx) ;
            numFeasPosGain[XDN] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(x1Ay1_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > threshold)
            {
                feasPosGains[XDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(x1Ay1_groups, &numIdx) ;
                numFeasPosGain[XDN] = numIdx ;
                if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
            }
        }
    }

    // y dnflip
    IdxGroupsByC_See_Best_Feasible_Cval(y1Ax2_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > threshold)
        {
            feasPosGains[YDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y1Ax2_groups, &numIdx) ;
            numFeasPosGain[YDN] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(y1Ax1_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > threshold)
            {
                feasPosGains[YDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(y1Ax1_groups, &numIdx) ;
                numFeasPosGain[YDN] = numIdx ;
                if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
            }
        }
    }

    VSPINT noFlipVal = strict ? 0 : -1 ;
    if(maxFeasPosGain == noFlipVal) // If maxFeasPosGain has its initial value
    {
        *retFlipVar = '\0' ;
        *retIndFlip = EMPTY ;
#ifdef DEBUG
        if(strict)
        {
            if(deep_search)
            {
                FlipSelector_Check_v5
                (
                    x, Ax, wgtx,
                    y, Ay, wgty,
                    *retFlipVar, *retIndFlip,
                    problem, IMPROVE_OBJ_ONLY_v5
                ) ;
            }
            else
            {
                FeasibleBeneficialFlipSelector_Check
                (
                    x, Ax, wgtx, y, Ay, wgty,
                    *retFlipVar, *retIndFlip, problem
                ) ;

                FeasibleStrictlyBeneficialFlipSelector_Check
                (
                    x, Ax, wgtx,
                    y, Ay, wgty,
                    *retFlipVar, *retIndFlip,
                    problem, IMPROVE_OBJ_ONLY
                ) ;
            }
        }
        else
        {
            FeasibleBeneficialFlipSelector_Check
            (
                x, Ax, wgtx, y, Ay, wgty,
                *retFlipVar, *retIndFlip, problem
            ) ;
        }
#endif
        return ;
    }

    /*------------------------------------------------------------------------*/
    // Record all flips with maximal gain in
    // maxGainFlips[4]     - a stack containing 0-3 which encode flips
    // m                   - size of maxGainFlips stack
    // numsMaxGainFlips[4] - a stack of the number of indices corr. to each flip
    // accum[5]            - running total of numsMaxGainFlips
    /*------------------------------------------------------------------------*/
    VSPINT maxGainFlips[4]     = {0} ;
    VSPINT m                   =  0  ;
    VSPINT numsMaxGainFlips[4] = {0} ;
    VSPINT acumm[5]            = {0} ;

    if( numFeasPosGain[XUP] > 0 && feasPosGains[XUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = XUP ;
        numsMaxGainFlips[m] = numFeasPosGain[XUP] ;
        m++ ;
        acumm[m] = acumm[m-1] + numFeasPosGain[XUP]  ;
    }


    if( numFeasPosGain[YUP] > 0 && feasPosGains[YUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = YUP ;
        numsMaxGainFlips[m] = numFeasPosGain[YUP] ;
        m++ ;
        acumm[m] = acumm[m-1] + numFeasPosGain[YUP]  ;
    }

    if( numFeasPosGain[XDN] > 0 && feasPosGains[XDN] == maxFeasPosGain )
    {
        maxGainFlips[m] = XDN ;
        numsMaxGainFlips[m] = numFeasPosGain[XDN] ;
        m++ ;
        acumm[m] = acumm[m-1] + numFeasPosGain[XDN]  ;
    }

    if( numFeasPosGain[YDN] > 0 && feasPosGains[YDN] == maxFeasPosGain )
    {
        maxGainFlips[m] = YDN ;
        numsMaxGainFlips[m] = numFeasPosGain[YDN] ;
        m++ ;
        acumm[m] = acumm[m-1] + numFeasPosGain[YDN] ;
    }

    /*------------------------------------------------------------------------*/
    // Choose a random flip with maximal gain
    // flipChoice - integer in [0, 3] representing the flip
    /*------------------------------------------------------------------------*/
    VSPINT baseRandNum = rand() ;
    VSPINT randChoice  = baseRandNum % acumm[m] ;
    VSPINT flipChoice ;
    bool was_initialized = false ;


    for(VSPINT p = 1 ; p <= m ; p++)
    {
        if(randChoice < acumm[p])
        {
            flipChoice = maxGainFlips[p-1] ;
            #ifdef DEBUG
                was_initialized = true ;
            #endif
            break ;
        }
    }

    #ifdef DEBUG
    if(!was_initialized)
    {
        ERR("Flip choice was not initialized!") ;
        fflush(stdout) ;
    }
    #endif

    #ifdef DEBUG
        if( !IN_RANGE(flipChoice, 0, 4) )
        {
            ERR("flipChoice incorrect") ;
        }
    #endif


    if(flipChoice == XUP)
    {
        IdxGroupsByC_See_Random_Idx_Given_C(x0Ay0_groups, baseRandNum, feasPosGains[X0], retIndFlip) ;
        *retFlipVar = 'x' ;
    }
    else if(flipChoice == YUP)
    {
        IdxGroupsByC_See_Random_Idx_Given_C(y0Ax0_groups, baseRandNum, feasPosGains[Y0], retIndFlip) ;
        *retFlipVar = 'y' ;
    }
    else if(flipChoice == XDN)
    {
        *retFlipVar = 'x' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C(x1Ay2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C(x1Ay1_groups, baseRandNum, retIndFlip) ;
        }
    }
    else
    {
        VSP_ASSERT(flipChoice == YDN) ;
        *retFlipVar = 'y' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C(y1Ax2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C(y1Ax1_groups, baseRandNum, retIndFlip) ;
        }
    }

    #ifdef DEBUG
        if(deep_search)
        {
            FlipSelector_Check_v5
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem,
                strict ? IMPROVE_OBJ_ONLY_v5 : IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5
            ) ;
        }
        else
        {
            FeasibleStrictlyBeneficialFlipSelector_Check
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem,
                strict ? IMPROVE_OBJ_ONLY : IMPROVE_OBJ_OR_NEUTRAL_DOWNFLIP
            ) ;
        }
    #endif

    return ;
}






// WARNING: This function is only for the strictly beneficial flips phase
// ASSUMED: The flip has NOT happened yet
// The "best feasible flip according to structs" is a flip which achieves the
// gain (accounting for feasibility)
// max{ c_(x0), c_(y0), gamma * p_(x1) - c_(x1), gamma * p_(y1) - c_(y1) }
// This function checks that the given flup achieves this gain
// If there is no flip with feasible gain according to the structs, then
// the function expects givenVarFlip = '\0' and givenIndFlip = EMPTY 
void FeasibleBeneficialFlipSelector_Check
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    char          givenVarFlip,
    VSPINT        givenIndFlip,
    Problem       *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, xi, yi, Axi, Ayi, ci ;
    VSPINT n = problem->graph->n ; 
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c = problem->graph->c ;
    VSPINT  u = problem->u ;


    /*------------------------------------------------------------------------*/
    // Get the p values
    /*------------------------------------------------------------------------*/
    VSPINT px0 = n + 1, py0 = n + 1, px1 = -1, py1 = -1 ;

    for(i = 0 ; i < n ; i++)
    {
        xi = x[i] ;
        yi = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;

        // Record px's
        if(xi)
        {
            if(Ayi > px1)
            {
                px1 = Ayi ;
            }
        }
        else
        {
            if(Ayi < px0)
            {
                px0 = Ayi ;
            }
        }

        // Record py's
        if(yi)
        {
            if(Axi > py1)
            {
                py1 = Axi ;
            }
        }
        else
        {
            if(Axi < py0)
            {
                py0 = Axi ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Record the c's and an index for each c
    // Also record:
    // x0witness, y0witness, x1witness, y1witness - indices achieving the
    // optimal p and c values
    /*------------------------------------------------------------------------*/
    VSPINT cx0 = -1, cy0 = -1, cx1 = gamma + 1, cy1 = gamma + 1 ;
    VSPINT x0witness = EMPTY, y0witness = EMPTY ;
    VSPINT x1witness = EMPTY, y1witness = EMPTY ;

    for(i = 0 ; i < n ; i++)
    {
        xi = x[i] ;
        yi = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;
        ci = c[i] ;

        // Deal with x
        if(xi == 0 && Ayi == px0)
        {
            if(ci > cx0 )
            {
                cx0 = ci ;
                x0witness = i ;
            }
        }

        if(xi == 1 && Ayi == px1)
        {
            if(ci < cx1 )
            {
                cx1 = ci ;
                x1witness = i ;
            }
        }

        // Deal with y
        if(yi == 0 && Axi == py0)
        {
            if(ci > cy0 )
            {
                cy0 = ci ;
                y0witness = i ;
            }
        }

        if(yi == 1 && Axi == py1)
        {
            if(ci < cy1 )
            {
                cy1 = ci ;
                y1witness = i ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Go through the 4 possible flips and record the maximum feasible gain
    /*------------------------------------------------------------------------*/
    VSPINT trueBestGain    = 0 ;
    VSPINT trueBestGainInd = EMPTY ;
    char   trueBestGainVar = '\0' ;
    VSPINT trueBestGainDir = 0 ; // 1 for upflip, -1 for downflip

    VSPINT x0gain = cx0 - gamma * px0 ;
    VSPINT y0gain = cy0 - gamma * py0 ;
    VSPINT x1gain = gamma * px1 - cx1 ;
    VSPINT y1gain = gamma * py1 - cy1 ;

    if( FEASIBLE_UPFLIP(wgtx, cx0, u) )
    {
        if(x0gain > trueBestGain)
        {
            trueBestGain    = x0gain ;
            trueBestGainInd = x0witness ;
            trueBestGainVar = 'x' ;
            trueBestGainDir = 1 ;
        }
    }

    if( FEASIBLE_UPFLIP(wgty, cy0, u) )
    {
        if(y0gain > trueBestGain)
        {
            trueBestGain    = y0gain ;
            trueBestGainInd = y0witness ;
            trueBestGainVar = 'y' ;
            trueBestGainDir = 1 ;
        }
    }


    if(x1gain > trueBestGain)
    {
        trueBestGain    = x1gain ;
        trueBestGainInd = x1witness ;
        trueBestGainVar = 'x' ;
        trueBestGainDir = -1 ;
    }

    if(y1gain > trueBestGain)
    {
        trueBestGain    = y1gain ;
        trueBestGainInd = y1witness ;
        trueBestGainVar = 'y' ;
        trueBestGainDir = -1 ;
    }


    /*------------------------------------------------------------------------*/
    // Compute the gain of the given flip. Load the given and best flips into
    // flip structs for easier debugging
    /*------------------------------------------------------------------------*/

    // Compute the gain of the given var
    VSPINT givenGain = EMPTY ;
    VSPINT dir ;

    if(givenIndFlip != EMPTY)
    {
        if(givenVarFlip == 'x')
        {
            dir = (1 - 2 * x[givenIndFlip]) ;
            givenGain = dir * (c[givenIndFlip] - gamma *(Ay[givenIndFlip] )) ;
        }
        else if(givenVarFlip == 'y')
        {
            dir = (1 - 2 * y[givenIndFlip]) ;
            givenGain = dir * (c[givenIndFlip] - gamma *(Ax[givenIndFlip] )) ;
        }
    }

    // Load the given flip into a Flip struct for easier debugging
    Flip givenFlip = {0} ;
    givenFlip.gamma = gamma ;
    givenFlip.i    = givenIndFlip ;
    givenFlip.var  = givenVarFlip ;
    givenFlip.gain = givenGain ;
    givenFlip.u = u ;
    if(givenVarFlip == 'x')
    {
        VSP_ASSERT(givenIndFlip != EMPTY) ;
        givenFlip.wgt = wgtx ;
        givenFlip.zi   = x[givenIndFlip] ;
        givenFlip.Awi  = Ay[givenIndFlip] ;
        givenFlip.ci   = c[givenIndFlip] ;
        if(x[givenIndFlip] == 0)
        {
            givenFlip.dir  = 1 ;
        }
        else
        {
            givenFlip.dir  = -1 ;
        }
    }
    else if(givenVarFlip == 'y')
    {
        VSP_ASSERT(givenIndFlip != EMPTY) ;
        givenFlip.wgt = wgty ;
        givenFlip.zi   = y[givenIndFlip] ;
        givenFlip.Awi  = Ax[givenIndFlip] ;
        givenFlip.dir  = 1 - 2 * y[givenIndFlip] ;
        givenFlip.ci   = c[givenIndFlip] ;
        if(y[givenIndFlip] == 0)
        {
            givenFlip.dir  = 1 ;
        }
        else
        {
            givenFlip.dir  = -1 ;
        }
    }

    // Load the true optimal flip into its own flip struct for easier debugging
    Flip trueFlip = {0} ;

    trueFlip.gamma = gamma ;
    trueFlip.i = trueBestGainInd ;
    trueFlip.gain = trueBestGain ;
    trueFlip.var = trueBestGainVar ;
    if(trueBestGainVar == 'x') trueFlip.wgt = wgtx ;
    else if(trueBestGainVar == 'y') trueFlip.wgt = wgty ;
    trueFlip.u = u ;

    if(trueBestGainInd != EMPTY)
    {
        trueFlip.ci = c[trueBestGainInd] ;
        
        if(trueBestGainVar == 'x')
        {
            trueFlip.zi  = x[trueBestGainInd] ;
            trueFlip.Awi = Ay[trueBestGainInd] ;
            trueFlip.dir = ( (1 - 2 * x[trueBestGainInd]) == 1 ) ? 1 : -1 ;
        }
        else
        {
            trueFlip.zi  = y[trueBestGainInd] ;
            trueFlip.Awi = Ax[trueBestGainInd] ;
            trueFlip.dir = ( (1 - 2 * y[trueBestGainInd]) == 1 ) ? 1 : -1 ;
        }
    }


    /*------------------------------------------------------------------------*/
    // Check the given flip against the true best flip
    /*------------------------------------------------------------------------*/
    // For debugging:
    // trueFlip struct
    // givenFlip struct

    VSPINT givenFlipExists  = (givenIndFlip != EMPTY) ;
    VSPINT trueBestExists   = (trueBestGainInd != EMPTY) ;

    // Either both must be true or both must be false simultaneously
    if( (givenFlipExists + trueBestExists != 0) && (givenFlipExists + trueBestExists != 2) ) 
    {
        ERR("Mismatch between given flip and true best flip") ;
    }

    if(trueBestExists)
    {
        if(givenFlip.gain != trueFlip.gain)
        {
            ERR("Mismatch between gain of given flip and true best flip") ;
        }
    }

    return ;
}



// Flip Selector Check




// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
// v2 differs from v1 only in how it breaks ties among max-gain flips:
// improving-imbalance flips are preferred over damaging-imbalance flips,
// and within each group downflips are preferred over upflips.
void SelectRandomOptimalBeneficialFlip_v2
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,
    bool deep_search,
    bool strict,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list, // only used in debugging

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list, // only used in debugging

    Problem       *problem,
    char           *retFlipVar, // '\0' if no flip
    VSPINT         *retIndFlip  // EMPTY if no flip
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        Check_All_Data
        (
            x, Ax, wgtx, y, Ay, wgty,
            x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
            y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
            problem
        ) ;
    #endif

    VSPINT u      = problem->u ;
    VSPINT gamma  = problem->graph->max_vtx_cost ;

    /*------------------------------------------------------------------------*/
    // Collect the following data for every feasible beneficial flip:
    // feasPosGains[k]   - gain of flip k (only valid if numFeasPosGain[k] > 0)
    // numFeasPosGain[k] - number of indices achieving that gain (0 if flip
    //                     k is not recorded; this doubles as a "present"
    //                     flag).
    // maxFeasPosGain    - max over recorded gains.
    //
    // A flip is "beneficial" / passes the threshold iff
    //   strict     : gain > 0
    //   non-strict : gain > -1   (i.e., gain >= 0)
    /*------------------------------------------------------------------------*/
    VSPINT cval, gain, numIdx ;
    VSPINT threshold      = strict ? 0 : -1 ;
    VSPINT maxFeasPosGain = threshold ;
    VSPINT feasPosGains  [4] = {0} ;
    VSPINT numFeasPosGain[4] = {0} ;

    // x upflip
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, deep_search, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > threshold)
        {
            feasPosGains[XUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(x0Ay0_groups, cval, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }

    // y upflip
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, deep_search, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > threshold)
        {
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(y0Ax0_groups, cval, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }

    // x downflip
    IdxGroupsByC_See_Best_Feasible_Cval(x1Ay2_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > threshold)
        {
            feasPosGains[XDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x1Ay2_groups, &numIdx) ;
            numFeasPosGain[XDN] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(x1Ay1_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > threshold)
            {
                feasPosGains[XDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(x1Ay1_groups, &numIdx) ;
                numFeasPosGain[XDN] = numIdx ;
                if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
            }
        }
    }

    // y dnflip
    IdxGroupsByC_See_Best_Feasible_Cval(y1Ax2_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > threshold)
        {
            feasPosGains[YDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y1Ax2_groups, &numIdx) ;
            numFeasPosGain[YDN] = numIdx ;
            if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(y1Ax1_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > threshold)
            {
                feasPosGains[YDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(y1Ax1_groups, &numIdx) ;
                numFeasPosGain[YDN] = numIdx ;
                if(gain > maxFeasPosGain) maxFeasPosGain = gain ;
            }
        }
    }

    VSPINT noFlipVal = strict ? 0 : -1 ;
    if(maxFeasPosGain == noFlipVal) // If maxFeasPosGain has its initial value
    {
        *retFlipVar = '\0' ;
        *retIndFlip = EMPTY ;
#ifdef DEBUG
        if(strict)
        {
            if(deep_search)
            {
                FlipSelector_Check_v5
                (
                    x, Ax, wgtx,
                    y, Ay, wgty,
                    *retFlipVar, *retIndFlip,
                    problem, IMPROVE_OBJ_ONLY_v5
                ) ;
            }
            else
            {
                FeasibleBeneficialFlipSelector_Check
                (
                    x, Ax, wgtx, y, Ay, wgty,
                    *retFlipVar, *retIndFlip, problem
                ) ;

                FeasibleStrictlyBeneficialFlipSelector_Check
                (
                    x, Ax, wgtx,
                    y, Ay, wgty,
                    *retFlipVar, *retIndFlip,
                    problem, IMPROVE_OBJ_ONLY
                ) ;
            }
        }
        else
        {
            FeasibleBeneficialFlipSelector_Check
            (
                x, Ax, wgtx, y, Ay, wgty,
                *retFlipVar, *retIndFlip, problem
            ) ;
        }
#endif
        return ;
    }

    /*------------------------------------------------------------------------*/
    // Partition the flips that achieve maxFeasPosGain into 4 priority tiers:
    //   tier 0: improving downflips  (reduces |wgtx - wgty|, is a downflip)
    //   tier 1: improving upflips    (reduces |wgtx - wgty|, is an upflip)
    //   tier 2: damaging  downflips  (fallback pool, downflips preferred)
    //   tier 3: damaging  upflips    (fallback pool)
    //
    // A flip "improves imbalance" iff
    //   - it is a downflip on the heavier side, or
    //   - it is an upflip   on the lighter side.
    // When wgtx == wgty no flip improves imbalance, so every max-gain flip
    // falls into tiers 2/3.
    //
    // The lowest-numbered non-empty tier is selected; within that tier, a
    // flip is chosen by weighted random selection (weighted by index count).
    //
    // tierFlips[t][.] - stack of flip codes in tier t
    // tierSizes[t]    - size of that stack
    // tierAcumm[t][.] - running totals of numFeasPosGain in that tier
    /*------------------------------------------------------------------------*/
    bool improves_imbalance[4] ;
    improves_imbalance[XUP] = (wgtx < wgty) ;
    improves_imbalance[YUP] = (wgty < wgtx) ;
    improves_imbalance[XDN] = (wgtx > wgty) ;
    improves_imbalance[YDN] = (wgty > wgtx) ;

    bool is_downflip[4] ;
    is_downflip[XUP] = false ;
    is_downflip[YUP] = false ;
    is_downflip[XDN] = true  ;
    is_downflip[YDN] = true  ;

    VSPINT tierFlips[4][4] = {{0}} ;
    VSPINT tierSizes[4]    = {0}   ;
    VSPINT tierAcumm[4][5] = {{0}} ;

    VSPINT allFlips[4] = {XUP, YUP, XDN, YDN} ;

    for(VSPINT i = 0 ; i < 4 ; i++)
    {
        VSPINT flip = allFlips[i] ;
        if(numFeasPosGain[flip] == 0)        continue ; // flip not recorded
        if(feasPosGains[flip] != maxFeasPosGain) continue ; // not at max

        VSPINT tier ;
        if(improves_imbalance[flip])
            tier = is_downflip[flip] ? 0 : 1 ;
        else
            tier = is_downflip[flip] ? 2 : 3 ;

        VSPINT slot = tierSizes[tier] ;
        tierFlips[tier][slot] = flip ;
        tierSizes[tier]++ ;
        tierAcumm[tier][slot + 1] = tierAcumm[tier][slot] + numFeasPosGain[flip] ;
    }

    // Pick the highest-priority (lowest-numbered) non-empty tier
    VSPINT chosenTier = EMPTY ;
    for(VSPINT t = 0 ; t < 4 ; t++)
    {
        if(tierSizes[t] > 0)
        {
            chosenTier = t ;
            break ;
        }
    }

    #ifdef DEBUG
        if(chosenTier == EMPTY)
        {
            ERR("No tier populated despite maxFeasPosGain being set!") ;
            fflush(stdout) ;
        }
    #endif

    /*------------------------------------------------------------------------*/
    // Choose a random flip within the chosen tier (weighted by index count)
    // flipChoice - integer in [0, 3] representing the flip
    /*------------------------------------------------------------------------*/
    VSPINT baseRandNum = rand() ;
    VSPINT randChoice  = baseRandNum % tierAcumm[chosenTier][tierSizes[chosenTier]] ;
    VSPINT flipChoice ;
    bool was_initialized = false ;

    for(VSPINT p = 1 ; p <= tierSizes[chosenTier] ; p++)
    {
        if(randChoice < tierAcumm[chosenTier][p])
        {
            flipChoice = tierFlips[chosenTier][p - 1] ;
            #ifdef DEBUG
                was_initialized = true ;
            #endif
            break ;
        }
    }

    #ifdef DEBUG
    if(!was_initialized)
    {
        ERR("Flip choice was not initialized!") ;
        fflush(stdout) ;
    }
    #endif

    #ifdef DEBUG
        if( !IN_RANGE(flipChoice, 0, 4) )
        {
            ERR("flipChoice incorrect") ;
        }
    #endif


    if(flipChoice == XUP)
    {
        IdxGroupsByC_See_Random_Idx_Given_C(x0Ay0_groups, baseRandNum, feasPosGains[X0], retIndFlip) ;
        *retFlipVar = 'x' ;
    }
    else if(flipChoice == YUP)
    {
        IdxGroupsByC_See_Random_Idx_Given_C(y0Ax0_groups, baseRandNum, feasPosGains[Y0], retIndFlip) ;
        *retFlipVar = 'y' ;
    }
    else if(flipChoice == XDN)
    {
        *retFlipVar = 'x' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C(x1Ay2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C(x1Ay1_groups, baseRandNum, retIndFlip) ;
        }
    }
    else
    {
        VSP_ASSERT(flipChoice == YDN) ;
        *retFlipVar = 'y' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C(y1Ax2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C(y1Ax1_groups, baseRandNum, retIndFlip) ;
        }
    }

    #ifdef DEBUG
        if(deep_search)
        {
            FlipSelector_Check_v5
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem,
                strict ? IMPROVE_OBJ_ONLY_v5 : IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5
            ) ;
        }
        else
        {
            FeasibleStrictlyBeneficialFlipSelector_Check
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem,
                strict ? IMPROVE_OBJ_ONLY : IMPROVE_OBJ_OR_NEUTRAL_DOWNFLIP
            ) ;
        }
    #endif

    return ;
}
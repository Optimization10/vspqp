#include "core.h"
#include "flips_private_v5.h"
#include "flips_common_api.h"


// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
void SelectRandomOptimalBeneficialFlip_v5
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,

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
    // feasPosGains[4] - in entry k, record the feasible positive gain (if any)
    //                  of the flip encoded as k
    // 
    // maxFeasPosGain - maximal feasible positive gain (weakly positive, >= 0)
    // numFeasPosGain - number of indices corresponding to that feasible
    //                  positive gain
    /*------------------------------------------------------------------------*/
    VSPINT cval, gain, numIdx ;
    VSPINT maxFeasPosGain = -1 ;
    VSPINT feasPosGains  [4] = {0} ;
    VSPINT numFeasPosGain[4] = {0} ;

    // x Upflip_v5
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, true, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain)
        {
            maxFeasPosGain = gain ;
            feasPosGains[XUP]   = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(x0Ay0_groups, cval, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
        }
    }

    // y Upflip_v5
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, true, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain)
        {
            maxFeasPosGain    = gain ;
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(y0Ax0_groups, cval, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
        }
    }

    // x Downflip_v5
    IdxGroupsByC_See_Best_Feasible_Cval(x1Ay2_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgtx, cval))
        {
            maxFeasPosGain = gain ;
            feasPosGains[XDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x1Ay2_groups, &numIdx) ;
            numFeasPosGain[XDN] = numIdx ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(x1Ay1_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgtx, cval))
            {
                maxFeasPosGain = gain ;
                feasPosGains[XDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(x1Ay1_groups, &numIdx) ;
                numFeasPosGain[XDN] = numIdx ;
            }
        }
    }

    // y dnflip
    IdxGroupsByC_See_Best_Feasible_Cval(y1Ax2_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgty, cval))
        {
            maxFeasPosGain = gain ;
            feasPosGains[YDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y1Ax2_groups, &numIdx) ;
            numFeasPosGain[YDN] = numIdx ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(y1Ax1_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgty, cval))
            {
                maxFeasPosGain = gain ;
                feasPosGains[YDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(y1Ax1_groups, &numIdx) ;
                numFeasPosGain[YDN] = numIdx ;
            }
        }
    }

    if(maxFeasPosGain == -1) // It's initialized as 0
    {
        *retFlipVar = '\0' ;
        *retIndFlip = EMPTY ;
        #ifdef DEBUG
            // Check_Given_Flip_Is_Optimal_v5
            // (
            //     x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
            // ) ;
            FlipSelector_Check_v5
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem, IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5
            ) ;
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

    if( feasPosGains[XUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = XUP ;
        numsMaxGainFlips[m] = numFeasPosGain[XUP] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[XUP]  ;
    }


    if( feasPosGains[YUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = YUP ;
        numsMaxGainFlips[m] = numFeasPosGain[YUP] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[YUP]  ;
    }    

    if( feasPosGains[XDN] == maxFeasPosGain )
    {
        maxGainFlips[m] = XDN ;
        numsMaxGainFlips[m] = numFeasPosGain[XDN] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[XDN]  ;
    }    

    if( feasPosGains[YDN] == maxFeasPosGain )
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


    for(VSPINT p = 1 ; p <= m ; p++)
    {
        if(randChoice < acumm[p])
        {
            flipChoice = maxGainFlips[p-1] ;
            break ;
        }
    }

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
        FlipSelector_Check_v5
        (
            x, Ax, wgtx,
            y, Ay, wgty,
            *retFlipVar, *retIndFlip,
            problem, IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5
        ) ;
    #endif

    return ;
}


// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
void SelectRandomOptimalStrictlyBeneficialFlip_v5
(
    VSPINT    *x,
    VSPINT    *y,
    VSPINT   *Ax,
    VSPINT   *Ay,
    VSPINT  wgtx,
    VSPINT  wgty,

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

    VSPINT u = problem->u ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    // Get a random number
    // VSPINT ranNum = rand() ;

    /*------------------------------------------------------------------------*/
    // Collect the following data
    // feasPosGains[4] - in entry k, record the feasible positive gain (if any)
    //                  of the flip encoded as k
    // 
    // maxFeasPosGain - maximal feasible positive gain
    // numFeasPosGain - number of indices corresponding to that feasible
    //                  positive gain
    /*------------------------------------------------------------------------*/
    VSPINT cval, gain, numIdx ;
    VSPINT maxFeasPosGain = 0 ;
    VSPINT feasPosGains  [4] = {0} ;
    VSPINT numFeasPosGain[4] = {0} ;

    // x Upflip_v5
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, true, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain)
        {
            maxFeasPosGain = gain ;
            feasPosGains[XUP]   = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(x0Ay0_groups, cval, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
        }
    }

    // y Upflip_v5
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, true, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain)
        {
            maxFeasPosGain    = gain ;
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(y0Ax0_groups, cval, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
        }
    }

    // x Downflip_v5
    // Do not have to do a deep search for downflips
    IdxGroupsByC_See_Best_Feasible_Cval(x1Ay2_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgtx, cval))
        {
            maxFeasPosGain = gain ;
            feasPosGains[XDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x1Ay2_groups, &numIdx) ;
            numFeasPosGain[XDN] = numIdx ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(x1Ay1_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgtx, cval))
            {
                maxFeasPosGain = gain ;
                feasPosGains[XDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(x1Ay1_groups, &numIdx) ;
                numFeasPosGain[XDN] = numIdx ;
            }
        }
    }

    // y dnflip
    IdxGroupsByC_See_Best_Feasible_Cval(y1Ax2_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain =  -(cval - gamma * 2) ;
        if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgty, cval))
        {
            maxFeasPosGain = gain ;
            feasPosGains[YDN] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y1Ax2_groups, &numIdx) ;
            numFeasPosGain[YDN] = numIdx ;
        }
    }
    else // if cval == EMPTY
    {
        IdxGroupsByC_See_Best_Feasible_Cval(y1Ax1_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
        if(cval != EMPTY)
        {
            gain =  -(cval - gamma * 1) ;
            if(gain > maxFeasPosGain && FEASIBLE_DNFLIP(wgty, cval))
            {
                maxFeasPosGain = gain ;
                feasPosGains[YDN] = gain ;
                IdxGroupsByC_Num_Idx_With_Top_C(y1Ax1_groups, &numIdx) ;
                numFeasPosGain[YDN] = numIdx ;
            }
        }
    }

    if(maxFeasPosGain == 0) // It's initialized as 0
    {
        *retFlipVar = '\0' ;
        *retIndFlip = EMPTY ;
        #ifdef DEBUG
            // Check_Given_Flip_Is_Optimal_v5
            // (
            //     x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
            // ) ;
            FlipSelector_Check_v5
            (
                x, Ax, wgtx,
                y, Ay, wgty,
                *retFlipVar, *retIndFlip,
                problem, IMPROVE_OBJ_ONLY
            ) ;
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

    if( feasPosGains[XUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = XUP ;
        numsMaxGainFlips[m] = numFeasPosGain[XUP] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[XUP]  ;
    }


    if( feasPosGains[YUP] == maxFeasPosGain )
    {
        maxGainFlips[m] = YUP ;
        numsMaxGainFlips[m] = numFeasPosGain[YUP] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[YUP]  ;
    }    

    if( feasPosGains[XDN] == maxFeasPosGain )
    {
        maxGainFlips[m] = XDN ;
        numsMaxGainFlips[m] = numFeasPosGain[XDN] ;
        m++ ; 
        acumm[m] = acumm[m-1] + numFeasPosGain[XDN]  ;
    }    

    if( feasPosGains[YDN] == maxFeasPosGain )
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
        // Check_Given_Flip_Is_Optimal_v5
        // (
        //     x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
        // ) ;

        FlipSelector_Check_v5
        (
            x, Ax, wgtx,
            y, Ay, wgty,
            *retFlipVar, *retIndFlip,
            problem, IMPROVE_OBJ_ONLY
        ) ;
        
    #endif

    return ;
}




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

        if(givenFlip.dir == -1 && !FEASIBLE_DNFLIP(givenFlip.wgt, givenFlip.ci) )
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
            if(gain > maxFeasGain && FEASIBLE_DNFLIP(wgtx, ci))
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
                trueFlip.dir = DOWN ;
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
            if(gain > maxFeasGain && FEASIBLE_DNFLIP(wgty, ci))
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
                trueFlip.dir = DOWN ;
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
    case IMPROVE_OBJ_ONLY :
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
    
    case GAIN_GEQ_0 :
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


void IdxGroupsByC_See_Top_Feasible_Cval_v5
(
    IdxGroupsByC      *z0Aw0_groups,
    Problem                *problem,
    VSPINT                     wgtz,
    VSPINT          *retMaxFeasCval,
    VSPINT                       *z, // Used for debugging only
    VSPINT                      *Aw // Used for debugging only
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        VSPINT n = problem->graph->n ;
        VSPINT maxFeasCval_dbg = EMPTY ;
        VSPINT witness = EMPTY ;
    #endif

    VSPINT  i, p, cval ;
    VSPINT  u = problem->u ;
    VSPINT *c = problem->graph->c ;

    VSPINT *heap = z0Aw0_groups->maxHeap->heap ;
    VSPINT sizeHeap = z0Aw0_groups->maxHeap->sizeHeap ;

    if(sizeHeap && FEASIBLE_UPFLIP(wgtz, heap[1], u) )
    {
        *retMaxFeasCval = heap[1] ;
        goto RET ;
    }

    VSPINT maxFeasCval = EMPTY ;
    for(p = 2 ; p <= sizeHeap ; p++)
    {
        if( FEASIBLE_UPFLIP(wgtz, heap[p], u) && heap[p] > maxFeasCval)
        {
            maxFeasCval = heap[p] ;
        }
    }
    *retMaxFeasCval = maxFeasCval ;
    
    RET:
    #ifdef DEBUG

        for(i = 0 ; i < n ; i++)
        {
            if( z[i] == 0 && Aw[i] == 0 && FEASIBLE_UPFLIP(wgtz, c[i], u) && c[i] > maxFeasCval_dbg)
            {
                maxFeasCval_dbg = c[i] ;
                witness = i ;
            }
        }

        if(*retMaxFeasCval != maxFeasCval_dbg)
        {
            ERR("Mismatch when searching for upflip") ;
        }
    #endif

    return ;
}
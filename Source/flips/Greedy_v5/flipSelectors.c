#include "core.h"
#include "flips_private_v5.h"
#include "flips_common_api.h"

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
    IdxGroupsByC_See_Top_Feasible_Cval_v5(x0Ay0_groups, problem, wgtx, &cval, x, Ay) ;
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
    IdxGroupsByC_See_Top_Feasible_Cval_v5(y0Ax0_groups, problem, wgty, &cval, y, Ax) ;
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
    IdxGroupsByC_See_Top_Cval_If_Feasible(x1Ay2_groups, wgtx, u, &cval) ;
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
        IdxGroupsByC_See_Top_Cval_If_Feasible(x1Ay1_groups, wgtx, u, &cval) ;
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
    IdxGroupsByC_See_Top_Cval_If_Feasible(y1Ax2_groups, wgty, u, &cval) ;
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
        IdxGroupsByC_See_Top_Cval_If_Feasible(y1Ax1_groups, wgty, u, &cval) ;
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
            IdxGroupsByC_See_Random_Idx_Top_C_v5(x1Ay2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(x1Ay1_groups, baseRandNum, retIndFlip) ;
        }
    }
    else
    {
        VSP_ASSERT(flipChoice == YDN) ;
        *retFlipVar = 'y' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(y1Ax2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(y1Ax1_groups, baseRandNum, retIndFlip) ;
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
    IdxGroupsByC_See_Top_Feasible_Cval_v5(x0Ay0_groups, problem, wgtx, &cval, x, Ay) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgtx, cval, u))
        {
            maxFeasPosGain = gain ;
            feasPosGains[XUP]   = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(x0Ay0_groups, cval, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
        }
    }

    // y Upflip_v5
    IdxGroupsByC_See_Top_Feasible_Cval_v5(y0Ax0_groups, problem, wgty, &cval, y, Ax) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgty, cval, u))
        {
            maxFeasPosGain    = gain ;
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Given_C(y0Ax0_groups, cval, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
        }
    }

    // x Downflip_v5
    IdxGroupsByC_See_Top_Cval_If_Feasible(x1Ay2_groups, wgtx, u, &cval) ;
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
        IdxGroupsByC_See_Top_Cval_If_Feasible(x1Ay1_groups, wgtx, u, &cval) ;
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
    IdxGroupsByC_See_Top_Cval_If_Feasible(y1Ax2_groups, wgty, u, &cval) ;
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
        IdxGroupsByC_See_Top_Cval_If_Feasible(y1Ax1_groups, wgty, u, &cval) ;
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
            IdxGroupsByC_See_Random_Idx_Top_C_v5(x1Ay2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(x1Ay1_groups, baseRandNum, retIndFlip) ;
        }
    }
    else
    {
        VSP_ASSERT(flipChoice == YDN) ;
        *retFlipVar = 'y' ;
        if(maxFeasPosGain >= gamma)
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(y1Ax2_groups, baseRandNum, retIndFlip) ;
        }
        else
        {
            IdxGroupsByC_See_Random_Idx_Top_C_v5(y1Ax1_groups, baseRandNum, retIndFlip) ;
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

void Compute_ps_And_cs_Dbg_v5
(
    VSPINT                    *x,
    VSPINT                   *Ax,
    VSPINT                  wgtx,
    VSPINT                    *y,
    VSPINT                   *Ay,
    VSPINT                  wgty,
    Problem         *problem,
    VSPINT                 ps[4],
    VSPINT                 cs[4],
    VSPINT          witnesses[4]
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, xi, yi, Axi, Ayi, ci ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;
    VSPINT  u    = problem->u ;

    /*------------------------------------------------------------------------*/
    // Get the p values
    /*------------------------------------------------------------------------*/
    VSPINT px0 = n + 1, py0 = n + 1, px1 = -1, py1 = -1 ;

    for(i = 0 ; i < n ; i++)
    {
        xi  = x[i] ;
        yi  = y[i] ;
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



    // Prepare return
    ps[X0] = px0 ;
    ps[X1] = px1 ;
    ps[Y0] = py0 ;
    ps[Y1] = py1 ;

    cs[X0] = cx0 ;
    cs[X1] = cx1 ;
    cs[Y0] = cy0 ;
    cs[Y1] = cy1 ;

    witnesses[X0] = x0witness ;
    witnesses[X1] = x1witness ;
    witnesses[Y0] = y0witness ;
    witnesses[Y1] = y1witness ;

    return ;
}

/*------------------------------------------------------------------------------
Encoding: 
0 --> x0
1 --> y0
2 --> x1
3 --> y1

inds[k] - k is for the encoding, inds[k] is a vertex number. inds[k] may be
          EMPTY

retIsFeas[k] - if inds[k] != EMPTY and the flip is feasible, the it's 1
               otherwise, it's 0

retGains[k] - if retIsFeas[k] = 1, then this is the gain of the flip. If
              retIsFeas[k] = 0, then this is garbage
------------------------------------------------------------------------------*/
void Compute_Feasibility_And_Gains_v5
(
    VSPINT                *Ax,
    VSPINT                *Ay,
    VSPINT               wgtx,
    VSPINT               wgty,
    Problem      *problem,
    VSPINT             inds[4],
    VSPINT        retIsFeas[4],
    VSPINT         retGains[4]
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i ;
    VSPINT n  = problem->graph->n ;
    VSPINT *c = problem->graph->c ;
    VSPINT u  = problem->u ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    for(i = 0 ; i < 4 ; i++) retIsFeas[i] = 0 ;

    if(inds[X0] != EMPTY && FEASIBLE_UPFLIP(wgtx, c[inds[X0]], u) )
    {
        retIsFeas[X0] = 1 ;
        retGains [X0] = GAIN_FLIP(1,  c[inds[X0]], gamma, Ay[inds[X0]]) ;
    }

    if(inds[Y0] != EMPTY && FEASIBLE_UPFLIP(wgty, c[inds[Y0]], u) )
    {
        retIsFeas[Y0] = 1 ;
        retGains [Y0] = GAIN_FLIP(1,  c[inds[Y0]], gamma, Ax[inds[Y0]]) ;
    }

    if(inds[X1] != EMPTY && FEASIBLE_DNFLIP(wgtx, c[inds[X1]]) )
    {
        retIsFeas[X1] = 1 ;
        retGains [X1] = GAIN_FLIP(-1,  c[inds[X1]], gamma, Ay[inds[X1]]) ;
    }

    if(inds[Y1] != EMPTY && FEASIBLE_DNFLIP(wgty, c[inds[Y1]]) )
    {
        retIsFeas[Y1] = 1 ;
        retGains [Y1] = GAIN_FLIP(-1,  c[inds[Y1]], gamma, Ax[inds[Y1]]) ;
    }

    return ;
}

void IdxGroupsByC_See_Top_Feasible_Cval_v5
(
    IdxGroupsByC      *z0Aw0_groups,
    Problem            *problem,
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
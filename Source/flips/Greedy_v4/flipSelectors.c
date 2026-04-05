#include "core.h"
#include "flips_private_v4.h"



// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
void SelectRandomOptimalBeneficialFlip
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

    // x upflip
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgtx, cval, u))
        {
            maxFeasPosGain = gain ;
            feasPosGains[XUP]   = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x0Ay0_groups, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
        }
    }

    // y upflip
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgty, cval, u))
        {
            maxFeasPosGain    = gain ;
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y0Ax0_groups, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
        }
    }

    // x downflip
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
            FeasibleBeneficialFlipSelector_Check
            (
                x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
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
    VSPINT flipChoice = NO_FLIP ;


    for(VSPINT p = 1 ; p <= m ; p++)
    {
        if(randChoice < acumm[p])
        {
            flipChoice = maxGainFlips[p-1] ;
            break ;
        }
    }

    #ifdef DEBUG
        if( !IN_RANGE(flipChoice, 0, 4) )
        {
            ERR("flipChoice incorrect") ;
        }
    #endif


    if(flipChoice == XUP)
    {
        IdxGroupsByC_See_Random_Idx_Top_C(x0Ay0_groups, baseRandNum, retIndFlip) ;
        *retFlipVar = 'x' ;
    }
    else if(flipChoice == YUP)
    {
        IdxGroupsByC_See_Random_Idx_Top_C(y0Ax0_groups, baseRandNum, retIndFlip) ;
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
        FeasibleStrictlyBeneficialFlipSelector_Check
        (
            x, Ax, wgtx,
            y, Ay, wgty,
            *retFlipVar, *retIndFlip,
            problem, IMPROVE_OBJ_OR_NEUTRAL_DOWNFLIP
        ) ;
    #endif

    return ;
}

// ASSUMED: p_(x, 1) and p_(y, 1) are each <= 2
void SelectRandomOptimalStrictlyBeneficialFlip
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

    // x upflip
    IdxGroupsByC_See_Best_Feasible_Cval(x0Ay0_groups, wgtx, u, false, &cval, x, Ay, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgtx, cval, u))
        {
            maxFeasPosGain = gain ;
            feasPosGains[XUP]   = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(x0Ay0_groups, &numIdx) ;
            numFeasPosGain[XUP] = numIdx ;
        }
    }

    // y upflip
    IdxGroupsByC_See_Best_Feasible_Cval(y0Ax0_groups, wgty, u, false, &cval, y, Ax, problem->graph->n, problem->graph->c) ;
    if(cval != EMPTY)
    {
        gain = cval ;
        if(gain > maxFeasPosGain && FEASIBLE_UPFLIP(wgty, cval, u))
        {
            maxFeasPosGain    = gain ;
            feasPosGains[YUP] = gain ;
            IdxGroupsByC_Num_Idx_With_Top_C(y0Ax0_groups, &numIdx) ;
            numFeasPosGain[YUP] = numIdx ;
        }
    }

    // x downflip
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
            FeasibleBeneficialFlipSelector_Check
            (
                x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
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
    VSPINT flipChoice = EMPTY ; ;


    for(VSPINT p = 1 ; p <= m ; p++)
    {
        if(randChoice < acumm[p])
        {
            flipChoice = maxGainFlips[p-1] ;
            break ;
        }
    }

    #ifdef DEBUG
        VSP_ASSERT(flipChoice != EMPTY) ;
    #endif

    if(flipChoice == XUP)
    {
        IdxGroupsByC_See_Random_Idx_Top_C(x0Ay0_groups, baseRandNum, retIndFlip) ;
        *retFlipVar = 'x' ;
    }
    else if(flipChoice == YUP)
    {
        IdxGroupsByC_See_Random_Idx_Top_C(y0Ax0_groups, baseRandNum, retIndFlip) ;
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
        FeasibleBeneficialFlipSelector_Check
        (
            x, Ax, wgtx, y, Ay, wgty, *retFlipVar, *retIndFlip, problem
        ) ;

        FeasibleStrictlyBeneficialFlipSelector_Check
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
void FeasibleStrictlyBeneficialFlipSelector_Check
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
    Expected_Flip_Effect expFlipEffect
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
    // Compute ps, cs, and witnesses (indices attaining the ps and cs)
    // For the witnesses, compute feasibility status, and gains
    /*------------------------------------------------------------------------*/

    VSPINT ps[4], cs[4], witness[4], isFeas[4], gains[4] ;
    Compute_ps_And_cs_Dbg
    (
        x, Ax, wgtx, NULL, y, Ay, wgty, NULL, problem,
        ps, cs, witness
    ) ;

    Compute_Feasibility_And_Gains
    (
        Ax, Ay, wgtx, wgty, problem,
        witness, isFeas, gains
    ) ;

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

    // In each case, trueFlip will be an optimal flip for that case. And we will
    // check givenFlip against it.
    // Initialize trueFlip as an empty flip
    Flip trueFlip = {0} ;
    trueFlip.var  = '\0' ;
    trueFlip.i    = EMPTY ;
    trueFlip.gamma = gamma ;
    trueFlip.u = u ;
    switch (expFlipEffect)
    {
    case IMPROVE_OBJ_ONLY :
        // Fill out trueFlip with a feasible, strictly beneficial flip if
        // possible
        bestGain = 0 ;
        for(k = 0 ; k < 4 ; k++)
        {
            if(isFeas[k] && gains[k] > bestGain)
            {
                i = witness[k] ;
                trueFlip.ci  = c[i] ;
                trueFlip.dir = GET_DIR(k) ;
                trueFlip.gain = gains[k] ;
                trueFlip.i = i ;
                trueFlip.zi = GET_ZI(k) ;
                if( GET_VAR(k) == X_VAR )
                {
                    trueFlip.var = 'x' ;
                    trueFlip.Awi = Ay[i] ;
                    trueFlip.wgt = wgtx ;
                }
                else
                {
                    trueFlip.var = 'y' ;
                    trueFlip.Awi = Ax[i] ;
                    trueFlip.wgt = wgty ;
                }
            }
        }

        if( !IFF(givenFlip.i != EMPTY, trueFlip.i != EMPTY) )
        {
            ERR("Mismatch in flips") ;
        }

        // Compare the given flip with true flip
        if(trueFlip.gain > 0)
        {
            if(givenFlip.var != '\0')
            {
                if(trueFlip.gain > givenFlip.gain)
                {
                    ERR("Suboptimal flip selected") ;
                }
            }
            else
            {
                ERR("Not detected a flip when there is one") ;
            }
        }
        else
        {
            if(givenFlip.i != EMPTY)
            {
                ERR("Seem to have detected a wrong flip") ;
            }
        }

        // Sanity check
        if(givenFlip.i != EMPTY && givenFlip.gain <= 0)
        {
            ERR("Bad flip selected, and code above is wrong") ;
        }

        break ;
    case IMPROVE_OBJ_OR_NEUTRAL_DOWNFLIP :
        // Get a true flip
        bestGain = 0 ;
        for(k = 0 ; k < 4 ; k++)
        {
            if(isFeas[k] && (gains[k] > bestGain || (k >= 2 && gains[k] >= bestGain) ) )
            {
                bestGain = gains[k] ;
                i = witness[k] ;
                trueFlip.ci  = c[i] ;
                trueFlip.dir = GET_DIR(k) ;
                trueFlip.gain = gains[k] ;
                trueFlip.i = i ;
                trueFlip.zi = GET_ZI(k) ;
                if( GET_VAR(k) == X_VAR )
                {
                    trueFlip.var = 'x' ;
                    trueFlip.Awi = Ay[i] ;
                    trueFlip.wgt = wgtx ;
                }
                else
                {
                    trueFlip.var = 'y' ;
                    trueFlip.Awi = Ax[i] ;
                    trueFlip.wgt = wgty ;
                }
            }
        }

        if(givenFlip.i != EMPTY)
        {
            // Check that it's either a strictly beneficial upflip or a neutral
            // downflip
            if(givenFlip.dir == 1)
            {
                if(givenFlip.gain <= 0)
                {
                    ERR("Flip does not strictly improve objective") ;
                }

                if(trueFlip.i == EMPTY || givenFlip.gain < trueFlip.gain)
                {
                    ERR("Error in flips") ;
                }
            }
            else
            {
                if(givenFlip.gain < 0)
                {
                    ERR("Flip damages objective") ;
                }
            }
        }
        else // no given flip
        {
            if(trueFlip.i != EMPTY)
            {
                ERR("Flip not found") ;
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
        // Get a true flip
        bestGain = -1 ;
        for(k = 0 ; k < 4 ; k++)
        {
            if(isFeas[k] && gains[k] > bestGain )
            {
                i = witness[k] ;
                trueFlip.ci  = c[i] ;
                trueFlip.dir = GET_DIR(k) ;
                trueFlip.gain = gains[k] ;
                trueFlip.i = i ;
                trueFlip.zi = GET_ZI(k) ;
                if( GET_VAR(k) == X_VAR )
                {
                    trueFlip.var = 'x' ;
                    trueFlip.Awi = Ay[i] ;
                    trueFlip.wgt = wgtx ;
                }
                else
                {
                    trueFlip.var = 'y' ;
                    trueFlip.Awi = Ax[i] ;
                    trueFlip.wgt = wgty ;
                }
            }
        }

        // Either both flips should exist or both should not exist
        if( !IFF(givenFlip.i != EMPTY, trueFlip.i != EMPTY) )
        {
            ERR("Mismatch between given and true flip") ;
        }

        // Compare the given flip with true flip
        if(trueFlip.i != EMPTY)
        {
            if(trueFlip.gain > givenFlip.gain)
            {
                ERR("Suboptimal flip selected") ;
            }
        }

        // Sanity check
        if(givenFlip.gain < 0)
        {
            ERR("Damaging flip selected, and code above is wrong") ;
        }
    }

    break ;
    default:
        ERR("Unknown expFlipEffect value") ;
        break;
    }
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

    if( FEASIBLE_DNFLIP(wgtx, cx1) )
    {
        if(x1gain > trueBestGain)
        {
            trueBestGain    = x1gain ;
            trueBestGainInd = x1witness ;
            trueBestGainVar = 'x' ;
            trueBestGainDir = -1 ;
        }
    }

    if( FEASIBLE_DNFLIP(wgty, cy1) )
    {
        if(y1gain > trueBestGain)
        {
            trueBestGain    = y1gain ;
            trueBestGainInd = y1witness ;
            trueBestGainVar = 'y' ;
            trueBestGainDir = -1 ;
        }
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
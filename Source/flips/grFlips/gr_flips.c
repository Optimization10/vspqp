#include "core.h"
#include "shared_api.h"
#include "math.h"
#include "score.h"
#include "states_manager.h"
#include "basic_flip_strats.h"
#include "flips_common_api.h"

// This code is finished except for calling flips
void GR_Flips(State *Xbest, State *Ybest, Problem *problem, VSPINT P)
{
    static long long unsigned count = 0 ;
    count++ ;

    IFDBG(check_states(Xbest, Ybest, problem)) ;
    /*------------------------------------------------------------------------*/
    // Local vars
    /*------------------------------------------------------------------------*/
    VSPINT p, k ;
    VSPFLOAT redFactor = 1.2;

    PointScore scoreBest, scoreWork ;

    bool success ;

    VSPFLOAT gammaRed,
    gammaOrg = problem->graph->max_vtx_cost ;

    StatesMan *SM    = StatesMan_Init(Xbest, Ybest) ;
    State *Xwork = StatesMan_GetState(SM, WORK_STATES, X_VAR) ;
    State *Ywork = StatesMan_GetState(SM, WORK_STATES, Y_VAR) ;
    IFDBG(check_states(Xbest, Ybest, problem)) ;
    IFDBG(check_states(Xwork, Ywork, problem)) ;

    FlipsCtx *FC = FlipsCtx_Init((State *[2]){Xwork, Ywork}, true, problem);

    /*------------------------------------------------------------------------*/
    // Main algorithm
    /*------------------------------------------------------------------------*/
    gammaRed = gammaOrg / redFactor ;
    success = true ;
    k = 0 ;
    PointScore_Update(&scoreBest, Xbest, Ybest, problem) ;
    for(p = 0 ; p < P ; p++)
    {
        IFDBG(check_states(Xbest, Ybest, problem)) ;
        IFDBG(check_states(Xwork, Ywork, problem)) ;
        if(!success)
        {
            StatesMan_UpdateWorkStates(SM) ;
            FlipsCtx_Set_Wgt(FC, 0, Xwork->weight) ;
            FlipsCtx_Set_Wgt(FC, 1, Ywork->weight) ;
            FlipsCtx_Initialize_Lists_And_Ts(FC) ;
        }

        // Apply feas ben flips to work states
        IFDBG(FlipsCtx_Check(FC)) ;
        BasicFlipStrats_FeasBeneficialFlips(FC, BOTH_VARS, false, gammaRed) ;

        // Apply feas ben flips to best states
        BasicFlipStrats_FeasBeneficialFlips(FC, BOTH_VARS, false, gammaOrg) ;

        // Compute workScore
        IFDBG(FlipsCtx_Check(FC)) ;
        Xwork->weight = FlipsCtx_Get_Wgt(FC, 0) ;
        Ywork->weight = FlipsCtx_Get_Wgt(FC, 1) ;
        PointScore_Update(&scoreWork, Xwork, Ywork, problem) ;
        success = PointScore_IsStrictlyBetter(&scoreWork, &scoreBest) ;

        if(success)
        {
            // printf("\nSuccess. Problem: %s  Level: %lld  n: %lld  objOld: %f  objNew: %f", problem->name, (long long) problem->level, (long long) problem->graph->n, scoreWork.obj, scoreBest.obj) ;
            StatesMan_UpdateBestStates(SM) ;
            gammaRed = gammaOrg / redFactor ;
            k = 0 ;
            scoreBest = scoreWork ;
        }
        else
        {
            // printf("\nFailure. Problem: %s  Level: %lld  n: %lld  objOld: %f  objNew: %f", problem->name, (long long) problem->level, (long long) problem->graph->n, scoreBest.obj, scoreWork.obj) ;
            gammaRed = gammaRed / redFactor ;
            k++ ;
        }
        IFDBG(check_states(Xbest, Ybest, problem)) ;
        IFDBG(check_states(Xwork, Ywork, problem)) ;
    }

    /*------------------------------------------------------------------------*/
    // Cleanup
    /*------------------------------------------------------------------------*/
    StatesMan_Destroy(SM) ;
    FlipsCtx_Destroy(FC) ;

    return ;
}
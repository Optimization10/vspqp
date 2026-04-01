#include "shared_api.h"
// #include "score.h"
// #include "breakpoints.h"
// #include "flips_context.h"
// #include "strategies.h"
// #include "states.h"

/*----------------------------------------------------------------------------*/
// FM Flips v8 module
// This is the main orchestrator for the FM flips v8 strategy
/*----------------------------------------------------------------------------*/
// void FM_Flips_v8
// (
//     VSP_state          *Zbest,
//     VSP_state          *Wbest,
//     VSP_problem      *problem,
//     VSPINT                  M  // Breakpoints at c^Tz = A_0 < A_1 ... < A_M = u.
//                                // If 0, then nothing happens 
// )
// {
//     SET_FUNCTION_COUNTER ;
//     VSPINT i, v, n = problem->Graph->n, bp, *lockedz ; 

//     // Initialize states manager, context, breakpoints
//     StatesMan   *SM = StatesMan_Init(Zbest, Wbest, problem) ;
//     FlipsCtx    *FC = FlipsCtx_Init((VSP_state*[]) { StatesMan_GetWork(SM, 0), StatesMan_GetWork(SM, 1) }, true, problem) ;
//     Breakpoints *BP = Breakpoints_Init(Zbest->weight, problem->u, M) ;

//     // Iterate over breakpoints
//     for(i = 1 ; i <= M ; i++)
//     {
//         // Lock all upflips in Z

//         // At this point, the point in flips context eneds to be a copy of
//         // the original point, FC needs to be up to data
//         // Lock all the downflips in Z

//         bp = Breakpoints_Get_Point(BP, i) ;

//         // Apply upflips until target weight is reached
//         Strategies_Do_Upflips_Until_Target_Wgt(FC, Z_VAR, true, FEAS_FLIP, bp) ;

//         // Apply all beneficial flips in W without locking
//         Strategies_Do_Beneficial_Flips(FC, W_VAR, false) ;

//         StatesMan_Update(SM) ;

//     }

// }
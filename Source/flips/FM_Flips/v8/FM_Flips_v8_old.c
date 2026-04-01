// #include "source/FM_Flips/V8/FM_Flips_v8.h"

// /*----------------------------------------------------------------------------*/
// /*----------------------------------------------------------------------------*/
// /*----------------------------------------------------------------------------*/
// // Let z be a variable NOT at the upper balance constraint (If both variables
// // are at the upper balance constraint of <= 0.5 and there are no constraint
// // violations, then we have a globally optimal solution. But this algorithms is
// // more for exploring situations where the graph is heavily unbalanced.

// // Assume c^Tz << u

// // Choose A_1, A_2, ..., A_M such that
// //     c^Tz A_0 = < A_1 < A_2 < ...< A_M < A_(M+1) = u
// // and these M + 2 points are equidistant.

// // Let d = (u - c^Tz) / M
// // For i = 0...M, A_i = c^Tz + i d

// // There should be two or three points floating around. A copy of the original,
// // the best point seen (which initially points to the original point), and a 
// // point in which we are making these changes.

// // Save a copy of the original point
// // For i = 1, 2, ..., m
// //     1) Starting from the original, lock entries i of z where z_i = 1
// //     2) Do upflips with locking in z until c^Tz just surpasses A_i
// //     3) Do all strictly beneficial flips in w
// //     4) Do all beneficial flips in (z, w)
// //     5) If the point reached is better than the original (in terms of
// //        objective or balance), then save it as the new best point
// // End For
// /*----------------------------------------------------------------------------*/
// /*----------------------------------------------------------------------------*/
// /*----------------------------------------------------------------------------*/


// // Zbest, Wbest - contains arrays of best point so far
// // XpermDeepCpy, YpermDeepCpy - a permanent deep copy of the input
// // Xwork, Ywork - At the beginning of each iteration, point is deeply copied
// //                from XpermDeepCpy
// // If at the end of the iteration, Xwork, Ywork is better than the best point,
// // then the contents of best states and work states are swapped

// void FM_Flips_v8
// (
//     VSP_state          *Zbest,
//     VSP_state          *Wbest,
//     VSP_problem      *problem,
//     VSPINT                  M  // Breakpoints at c^Tz = A_0 < A_1 ... < A_M = u.
//                                // If 0, then nothing happens 
// )
// {
//     if(M == 0) return ;
//     /*------------------------------------------------------------------------*/
//     // Debugging checks for input
//     /*------------------------------------------------------------------------*/
//     static long long unsigned count = 0 ;
//     count++ ;
//     #ifdef VSP_DEBUG
//         // Check input data
//         VSP_check_states(Zbest, Wbest, problem) ;

//         // Check that x and y are binary
//         if(Zbest->point->ind_frac != EMPTY || Wbest->point->ind_frac != EMPTY)
//         {
//             ERR("Non-binary points in FM flips") ;
//         }
        
//         // Record initial objective
//         VSPINT objInitDbg ;
//         Compute_ObjValSimple
//         (
//             Zbest->point->int_part, Zbest->Az->int_part, Zbest->weight,
//             Wbest->point->int_part, Wbest->Az->int_part, Wbest->weight,
//             problem, &objInitDbg
//         ) ;
//     #endif


//     /*------------------------------------------------------------------------*/
//     // Create states for juggling
//     /*------------------------------------------------------------------------*/
//     VSPINT v, objBest, imbBest, i, d, Ai, objWork, imbWork, zi,
//     n = problem->Graph->n,
//     u = problem->u ;

//     // Best point: compute objective and imbalance
//     VSP_state *statesBest[2] = {Zbest, Wbest} ;
//     Compute_ObjValSimple
//     (
//         statesBest[Z_VAR]->point->int_part, statesBest[Z_VAR]->Az->int_part, statesBest[Z_VAR]->weight,
//         statesBest[W_VAR]->point->int_part, statesBest[W_VAR]->Az->int_part, statesBest[W_VAR]->weight,
//         problem, &objBest
//     ) ;
//     imbBest = IMBALANCE(Zbest->weight, Wbest->weight) ;

//     // Permanent copy of original point: allocate and copy
//     VSP_state *statesOrgCpy[2] ;
//     for(v = 0 ; v < 2 ; v++) statesOrgCpy[v] = VSP_construct_state(n) ;
//     VSP_copy_state( Zbest, statesOrgCpy[Z_VAR], problem ) ;
//     VSP_copy_state( Wbest, statesOrgCpy[W_VAR], problem ) ;

//     // Work states: allocate
//     VSP_state *statesWork[2];
//     for(v = 0 ; v < 2 ; v++) statesWork[v] = VSP_construct_state(n) ;


//     /*------------------------------------------------------------------------*/
//     // Create As
//     /*------------------------------------------------------------------------*/



//     /*------------------------------------------------------------------------*/
//     // Create linked lists, Ts, and locked arrays
//     /*------------------------------------------------------------------------*/
//     VSPINT k, indFlip, varFlip, gainFlip,
//     nflips = 0,
//     *locked[2] ;


//     // Linked Lists
//     Doubly_LLs_Sized *Lists[2][2] ;
//     for(v = 0 ; v < 2 ; v++)
//     {
//         for(i = 0 ; i < 2 ; i++)
//         {
//             Lists[v][i] = Doubly_LLs_Sized_Init(n) ;
//         }
//     }

//     VSPINT Ts[2][2] ;

//     // Locked arrays
//     for(v = 0 ; v < 2 ; v++)
//     {
//         locked[v] = malloc(n * sizeof(VSPINT)) ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Do Flips
//     /*------------------------------------------------------------------------*/
//     for(i = 1 ; i <= M ; i++)
//     {        
//         // Deep copy the original point into the work point
//         VSP_copy_state(statesOrgCpy[Z_VAR], statesWork[Z_VAR], problem) ;
//         VSP_copy_state(statesOrgCpy[W_VAR], statesWork[W_VAR], problem) ;

//         /*--------------------------------------------------------------------*/
//         // Place data into 0-1 format
//         /*--------------------------------------------------------------------*/
//         VSPINT
//         *vars [2] = { statesWork[Z_VAR]->point->int_part  ,   statesWork[W_VAR]->point->int_part  },
//         *Avars[2] = { statesWork[Z_VAR]->Az->int_part     ,   statesWork[W_VAR]->Az->int_part     },
//         *wgts [2] = { &(statesWork[Z_VAR]->weight)        ,   &(statesWork[W_VAR]->weight)        }  ;

//         // Randomly initialize data structures
//         for(v = 0 ; v < 2 ; v++)
//         {
//             Locked_Arr_Reset(locked[v], n) ;
//         }
//         Initialize_Lists_And_Ts
//         (
//             vars, Avars, locked, Lists, Ts, problem, TRUE
//         ) ;
//         IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;


//         // Lock entries of i with z_i = 1
//         for(k = 0 ; k < n ; k++)
//         {
//             if(vars[Z_VAR][k] == 1) locked[Z_VAR][k]++ ;
//             // Doubly_LLs_Sized_Delete(Lists[Z_VAR], )
//         }        

//         Ai = As[i] ;
//         Do_Upflips_Until_Target_Wgt
//         (
//             vars, Avars, wgts, locked, Ts, Lists, problem,
//             Z_VAR, true, FEAS_FLIP, Ai
//         ) ;

//         // Do all beneficial flips in w, without locking
//         // (so we don't have to unlock later)
//         Do_Beneficial_Flips
//         (
//             vars, Avars, wgts, locked, Ts, Lists, problem,
//             W_VAR, false
//         ) ;

//         // Do all beneficial flips in (z, w)
//         Do_Beneficial_Flips
//         (
//             vars, Avars, wgts, locked, Ts, Lists, problem,
//             BOTH_VARS, false
//         ) ;

//         // If the point reached is better than the original, save it as the new
//         // best point by swapping the contents of the work point with the best
//         // point
//         Compute_ObjValSimple
//         (
//             statesWork[Z_VAR]->point->int_part, statesWork[Z_VAR]->Az->int_part, statesWork[Z_VAR]->weight,
//             statesWork[W_VAR]->point->int_part, statesWork[W_VAR]->Az->int_part, statesWork[W_VAR]->weight,
//             problem, &objWork
//         ) ;
//         imbWork = IMBALANCE(statesWork[Z_VAR]->weight, statesWork[W_VAR]->weight) ;

//         if(objWork > objBest || (objWork == objBest && imbWork < imbBest) )
//         {
//             // Swap work states with best states
//             for(v = 0 ; v < 2 ; v++) Swap_Mem(statesBest[v], statesWork[v], sizeof(VSP_state)) ;
//             objBest = objWork ;
//             imbBest = imbWork ; 
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Free memory
//     /*------------------------------------------------------------------------*/
//     // Free LLs
//     for(v = 0 ; v < 2 ; v++)
//     {
//         for(i = 0 ; i < 2 ; i++) Doubly_LLs_Sized_Destroy(Lists[v][i]) ;
//     }


//     // Free juggling states
//     for(v = 0 ; v < 2 ; v++)
//     {
//         VSP_destruct_state( statesOrgCpy[v] ) ;
//         VSP_destruct_state( statesWork[v] ) ;
//     }

//     // Free locked arrays
//     for(v = 0 ; v < 2 ; v++) free(locked[v]) ;
    

//     free(As) ;

//     /*------------------------------------------------------------------------*/
//     // Debugging checks for output
//     /*------------------------------------------------------------------------*/

//     return ;
// }









// static void Locked_Arr_Reset(VSPINT *locked, VSPINT n)
// {
//     SET_FUNCTION_COUNTER ;
//     VSP_ASSERT(locked != NULL) ;
//     for(VSPINT i = 0 ; i < n ; i++) locked[i] = 0 ;
//     return ;
// }
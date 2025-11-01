#include "shared_api.h"
// #include "source/FM_Flips/FM_Flips_v8.h"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Let z be a variable NOT at the upper balance constraint (If both variables
// are at the upper balance constraint of <= 0.5 and there are no constraint
// violations, then we have a globally optimal solution. But this algorithms is
// more for exploring situations where the graph is heavily unbalanced.

// Assume c^Tz << u

// Choose A_1, A_2, ..., A_M such that
//     c^Tz A_0 = < A_1 < A_2 < ...< A_M < A_(M+1) = u
// and these M + 2 points are equidistant.

// Let d = (u - c^Tz) / M
// For i = 0...M, A_i = c^Tz + i d

// There should be two or three points floating around. A copy of the original,
// the best point seen (which initially points to the original point), and a 
// point in which we are making these changes.

// Save a copy of the original point
// For i = 1, 2, ..., m
//     1) Starting from the original, lock entries i of z where z_i = 1
//     2) Do upflips with locking in z until c^Tz just surpasses A_i
//     3) Do all strictly beneficial flips in w
//     4) Do all beneficial flips in (z, w)
//     5) If the point reached is better than the original (in terms of
//        objective or balance), then save it as the new best point
// End For
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/


// Zbest, Wbest - contains arrays of best point so far
// XpermDeepCpy, YpermDeepCpy - a permanent deep copy of the input
// Xwork, Ywork - At the beginning of each iteration, point is deeply copied
//                from XpermDeepCpy
// If at the end of the iteration, Xwork, Ywork is better than the best point,
// then the contents of best states and work states are swapped

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
//     VSP_copy_state( Zbest, statesOrgCpy[Z_VAR] ) ;
//     VSP_copy_state( Wbest, statesOrgCpy[W_VAR] ) ;

//     // Work states: allocate
//     VSP_state *statesWork[2];
//     for(v = 0 ; v < 2 ; v++) statesWork[v] = VSP_construct_state(n) ;


//     /*------------------------------------------------------------------------*/
//     // Create As
//     /*------------------------------------------------------------------------*/
//     VSPINT *As = malloc( (M + 1) * sizeof(VSPINT)) ;
//     d = (u - Zbest->weight) / CAST(M, VSPFLOAT) ;
//     for(i = 0 ; i < M ; i++)
//     {
//         As[i] = Zbest->weight + d * i ; // rounds down
//         VSP_ASSERT(As[i] <= u) ;
//     }
//     As[M] = u ;


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
//         VSP_copy_state(statesOrgCpy[Z_VAR], statesWork[Z_VAR]) ;
//         VSP_copy_state(statesOrgCpy[W_VAR], statesWork[W_VAR]) ;

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


// // whichVars = 0 --- var 0
// // whichVars = 1 --- var 1
// // whichVars = 2 --- both vars

// // mode = 0 --- select the best feasible, even with negative gain 
// // mode = 1 --- select the best beneficial flip. Return EMPTY if no
// //              feas ben flip
// static void Select_Flip
// (
//     VSPINT                   *vars[2],
//     VSPINT                  *Avars[2],
//     VSPINT                   *wgts[2],
//     VSPINT                 *locked[2],
//     VSPINT                   Ts[2][2],
//     Doubly_LLs_Sized     *Lists[2][2],
//     VSP_problem              *problem,
//     VSPINT                *retFlipVar,
//     VSPINT                *retIndFlip,
//     VSPINT               *retFlipGain,
//     Flip_Selector_Vars      whichVars,
//     Flip_Selector_Mode           mode
// )
// {
//     SET_FUNCTION_COUNTER ;
//     IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
//     VSP_ASSERT( IN_RANGE(whichVars, 0, 3) ) ;

//     VSPINT  i, v, zi, k, ci, T, feasGains[4], start, end,
//     flipInds [4] = {EMPTY, EMPTY, EMPTY, EMPTY},
//     numFlips [4] = {0},
//     gamma = problem->Graph->max_vtx_cost,
//     n     = problem->Graph->n,
//     u     = problem->u,
//     *c    = problem->Graph->c ;

//     /*------------------------------------------------------------------------*/
//     // Collect the following data from the variables specified by whichVars
//     // flipInds[4]  - in entry k, record the index of the best flip of type k
//     //                  if such a flip exists and is feasible. Otherwise,
//     //                  leave the initial value of -1 unchanged
//     // feasGains[4] - in entry k, record the feasible gain (if any)
//     //                  of the flip encoded as k
//     // numFlips[4]  - in entry k, if flip k exists and is feasible, store 
//     //                the size of the optimal heap
//     // maxFeasGain  - maximal feasible positive gain (weakly positive, >= 0)
//     /*------------------------------------------------------------------------*/
//     const VSPINT noGainFlag = -gamma * n - gamma  ;
//     VSPINT maxFeasGain = noGainFlag ;

//     // 0 --> (0, 1)
//     // 1 --> (1, 2)
//     // 2 --> (0, 2)
//     start = whichVars % 2 ;
//     end   = -(whichVars * whichVars - 3 * whichVars) / 2 + 1 ;
//     for(v = start ; v < end ; v++)
//     {
//         for(zi = 0 ; zi < 2 ; zi++)
//         {
//             T = Ts[v][zi] ;
//             // If looking for feasible flips, skip upflips T > 0
//             if(mode == FEAS_BEN_FLIP && zi == 0 && T != 0) continue ;

//             if( IS_VALID_T(zi, T) )
//             {
//                 Doubly_LLs_Sized_See_Entry(Lists[v][zi], T, &i) ;
//                 VSP_ASSERT(i != EMPTY) ;

//                 // If flip is feasible, compute gain and store it in feasGains
//                 ci = c[i] ;
//                 if( IS_FEASIBLE_FLIP(zi, *wgts[v], ci, u ) )
//                 {
//                     k = GET_ENCODING(v, zi) ;
//                     flipInds [k] = i ;
//                     numFlips [k] = Lists[v][zi]->sizes[T] ;
//                     feasGains[k] = GAIN_FLIP( DIR_FLIP(zi), ci, gamma, T ) ;
//                     if( feasGains[k] > maxFeasGain )
//                     {
//                         maxFeasGain = feasGains[k] ;
//                     } 
//                 }
//             }
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Quit early if there is no feasible gain
//     /*------------------------------------------------------------------------*/
//     if(maxFeasGain == noGainFlag)
//     {
//         *retFlipVar  = NO_VAR ;
//         *retIndFlip  = EMPTY ;
//         *retFlipGain = maxFeasGain ;
//         IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
//         return ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Record all flips with maximal gain in 
//     // maxGainFlips[4]     - a stack containing 0-3 which encode flips
//     // m                   - size of maxGainFlips stack
//     // numsMaxGainFlips[4] - a stack of the number of indices corr. to each flip
//     // accum[5]            - running total of numsMaxGainFlips
//     /*------------------------------------------------------------------------*/
//     VSPINT maxGainFlips[4]     = {0} ;
//     VSPINT m                   =  0  ;
//     VSPINT numsMaxGainFlips[4] = {0} ;
//     VSPINT acumm[5]            = {0} ;

//     for(k = 0 ; k < 4 ; k++)
//     {
//         if(flipInds[k] != EMPTY && feasGains[k] == maxFeasGain)
//         {
//             maxGainFlips    [m] = k ;
//             numsMaxGainFlips[m] = numFlips[k] ;
//             m++ ; 
//             acumm[m]            = acumm[m-1] + numFlips[k] ;
//         }
//     }
//     VSP_ASSERT(m >= 1) ;


//     /*------------------------------------------------------------------------*/
//     // Choose a random flip with maximal gain
//     // flipChoice - integer in [0, 3] representing the flip
//     /*------------------------------------------------------------------------*/
//     VSPINT randChoice  = rand() % acumm[m] ;
//     VSPINT flipChoice ;


//     for(VSPINT p = 1 ; p <= m ; p++)
//     {
//         if(randChoice < acumm[p])
//         {
//             flipChoice = maxGainFlips[p-1] ;
//             break ;
//         }
//     }
//     *retFlipVar  = GET_VAR(flipChoice) ;
//     *retIndFlip  = flipInds[flipChoice] ;
//     *retFlipGain = maxFeasGain ;

//     VSPINT dir = GET_DIR(flipChoice) ;

//     IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
//     return ;
// }




// // Let z be the variable of the flip, and w the other variable
// // ASSUMED: The flip is happening with an Awi or Azi with optimal T
// static void Do_Flip
// (
//     VSPINT               *vars[2],
//     VSPINT              *Avars[2],
//     VSPINT               *wgts[2],
//     VSPINT             *locked[2],
//     VSPINT               Ts[2][2],
//     Doubly_LLs_Sized *Lists[2][2],
//     VSPINT                    var, // variable of the flip
//     VSPINT                      i, // index of the flip
//     VSP_problem          *problem,
//     bool                   lockIt  // lock the flip we do
// )
// {
//     SET_FUNCTION_COUNTER ;
//     IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

//     VSPINT p, j, T, wj, Azj, ziNew,
//     otherVar = OTHER_VAR(var) ;

//     VSPINT
//     n   = problem->Graph->n,
//     *z  = vars [var],
//     *Az = Avars[var],
//     *w  = vars [otherVar],
//     *Aw = Avars[otherVar],
//     zi  = vars [var][i] ;

//     VSPINT
//     dir = DIR_FLIP(zi),
//     Awi = Ts[var][zi],
//     *c  = problem->Graph->c,
//     *Ai = problem->Graph->Ai,
//     *Ap = problem->Graph->Ap ;

//     Doubly_LLs_Sized *lists ;

//     /*------------------------------------------------------------------------*/
//     // Remove i from its list
//     // Udate the T corresponding to zi
//     // Update weights
//     // Update zi to di + dir and lock zi
//     /*------------------------------------------------------------------------*/
//     Doubly_LLs_Sized_Delete(Lists[var][zi], Awi, i) ;

//     // If i was the only entry in its list, then the T needs to be updated
//     T = Awi ;
//     if( zi ) // zi = 1
//     {
//         // We are interested in max, so search downwards
//         while(T > -1 && Lists[var][zi]->sizes[T] == 0 ) T-- ;
//     }
//     else
//     {
//         // We are interested in min, so search upwards
//         while(T < n + 1 && Lists[var][zi]->sizes[T] == 0 ) T++ ;
//     }

//     // Update T
//     Ts[var][zi] = T ;

//     // Update zi, weight, and locked status
//     z[i]       += dir ;
//     *wgts[var] += dir * c[i] ;

//     if(lockIt) locked[var][i]++ ;

//     /*------------------------------------------------------------------------*/
//     // If not locking zi:
//     // Add i to new list
//     // Udate the T corresponding to new zi
//     /*------------------------------------------------------------------------*/
//     ziNew = zi + dir ;
//     if(!lockIt)
//     {
//         Doubly_LLs_Sized_Insert(Lists[var][ ziNew ], Awi, i) ;

//         // If Awi is better than the current T
//         T = Ts[var][ziNew] ;
//         if(dir * Awi > dir * T) Ts[var][ziNew]  = Awi ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Loop through j in N(i) cup (i) where wj is unlocked, and shift them one
//     // list up or down
//     /*------------------------------------------------------------------------*/
//     VSPLONG start = Ap[i] ;
//     VSPLONG end   = Ap[i+1] ;
//     p = start ;
//     j = i ;

//     bool iItselfFinished = false ;
//     do
//     {
//         if(iItselfFinished)
//         {
//             j = Ai[p] ;
//             p++ ;
//         }
//         iItselfFinished = true ;

//         // if wj is locked, do nothing
//         if( !locked[otherVar][j] )
//         {
//             // Remove j from its list
//             wj    = w[j] ;
//             Azj   = Az[j] ;
//             lists = Lists[otherVar][wj] ;

//             Doubly_LLs_Sized_Delete(lists, Azj, j) ;
//             Doubly_LLs_Sized_Insert(lists, Azj + dir, j) ;

//             // If wj = 0 and dir = 1, or if wj = 1 and dir = -1
//             if(dir + 2 * wj == 1)
//             {
//                 T = Ts[otherVar][wj] ;
//                 if(lists->sizes[T] == 0 )
//                 {
//                     Ts[otherVar][wj] += dir ;
//                 }
//             }
//             else // If wj = 0 and dir = -1, or if wj = 1 and dir = 1
//             {
//                 if( Azj == Ts[otherVar][wj])
//                 {
//                     Ts[otherVar][wj] += dir ;
//                 }
//             }     
//         }
//         Az[j] += dir ;
//     } while ( p != end ) ;

//     IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
//     return ;
// }



// static void Locked_Arr_Reset(VSPINT *locked, VSPINT n)
// {
//     SET_FUNCTION_COUNTER ;
//     VSP_ASSERT(locked != NULL) ;
//     for(VSPINT i = 0 ; i < n ; i++) locked[i] = 0 ;
//     return ;
// }


// // Initialize sized doubly linked lists and Ts. Ignore locked vertices. If
// // randomize is true, then build the lists in a randomized manner by going
// // through the vertices in a random order. Otherwise, go through the vertices
// // in numerical order
// // locked - If null, then ignored. One or more of its entries can also be NULL,
// //          in which case that particular array is ignored
// void Initialize_Lists_And_Ts
// (
//     VSPINT                 *vars[2],
//     VSPINT                *Avars[2],
//     VSPINT               *locked[2],
//     Doubly_LLs_Sized   *lists[2][2],
//     VSPINT                 Ts[2][2],
//     VSP_problem            *problem,
//     VSPINT                randomize  // 0 or 1
// )
// {
//     SET_FUNCTION_COUNTER ;
//     VSPINT i, v, j,*perm, zi, Awi, encoding, dir,
//     n = problem->Graph->n ;

//     /*------------------------------------------------------------------------*/
//     // Reset lists
//     /*------------------------------------------------------------------------*/
//     for(v = 0 ; v < 2 ; v++)
//     {
//         for(zi = 0 ; zi < 2 ; zi++)
//         {
//             Doubly_LLs_Sized_Reset(lists[v][zi]) ;
//         }
//     }


//     if(randomize)
//     {
//         perm = malloc(n * sizeof(*perm)) ;
//         Get_Permutation(perm, n) ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Initialize
//     /*------------------------------------------------------------------------*/
//     // TODO: If no heap, use D+2 where D is max degree
//     for(v = 0 ; v < 2 ; v++)
//     {
//         Ts[v][0] = n + 1 ;
//         Ts[v][1] = -1 ;
//     }

//     for(j = 0 ; j < n ; j++)
//     {
//         i = randomize ? perm[j] : j ;
//         for(v = 0 ; v < 2 ; v++)
//         {
//             // Skip locked vertices
//             if( !UNLOCKED(locked[v], i) ) continue ;

//             zi  = vars[v][i] ;
//             Awi = Avars[GET_OTHER_VAR(v)][i] ;
//             dir = DIR_FLIP(zi) ;

//             // Insert into ll and update T
//             Doubly_LLs_Sized_Insert(lists[v][zi], Awi, i);
//             if(dir * Awi < dir * Ts[v][zi] ) Ts[v][zi] = Awi ;
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Cleanup
//     /*------------------------------------------------------------------------*/
//     if(randomize)
//     {
//         free(perm) ;
//     }
//     return ;
// }



// static void Check_Data
// (
//     VSPINT               *vars[2],
//     VSPINT              *Avars[2],
//     VSPINT               *wgts[2],
//     VSPINT             *locked[2],
//     VSPINT               Ts[2][2],
//     Doubly_LLs_Sized *Lists[2][2],
//     VSP_problem          *problem
// )
// {
//     static long long unsigned count = 0 ;
//     count++ ;

//     WARN_IF_NOT_DEBUG_MODE ;

//     VSPINT v, zi, T, Tdbg, k,
//     n = problem->Graph->n ;
//     /*------------------------------------------------------------------------*/
//     // Check vars, Avars, weights, and feasibility by placing into states
//     /*------------------------------------------------------------------------*/
//     VSP_point
//     ptX = {.int_part = vars[Z_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n},
//     ptY = {.int_part = vars[W_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n} ;

//     VSP_Az
//     structAx = {.int_part   = Avars[Z_VAR] , .ind_frac_z = EMPTY , .n = n   ,
//                 .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   },
//     structAy = {.int_part   = Avars[W_VAR] , .ind_frac_z = EMPTY , .n = n   ,
//                 .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   };

//     VSP_state
//     X = {.point = &ptX, .Az = &structAx, .weight = *wgts[Z_VAR]},
//     Y = {.point = &ptY, .Az = &structAy, .weight = *wgts[W_VAR]} ;
//     VSP_check_states(&X, &Y, problem) ;

//     /*------------------------------------------------------------------------*/
//     // Check linked lists
//     /*------------------------------------------------------------------------*/
//     for(v = 0 ; v < 2 ; v++)
//     {
//         Doubly_LLs_Sized_Check
//         (
//             Lists[v][0], Lists[v][1] , vars[v], Avars[OTHER_VAR(v)], locked[v]
//         ) ;
//     }

//     /*------------------------------------------------------------------------*/
//     // Check Ts
//     /*------------------------------------------------------------------------*/
//     VSPINT ps[4], cs[4], wit[4] ;
//     Compute_ps_And_cs_Dbg
//     (
//         vars[Z_VAR], Avars[Z_VAR], *wgts[Z_VAR], locked[Z_VAR],
//         vars[W_VAR], Avars[W_VAR], *wgts[W_VAR], locked[W_VAR],
//         problem, ps, cs, wit
//     ) ;

//     for(v = 0 ; v < 2 ; v++)
//     {
//         for(zi = 0 ; zi < 2 ; zi++)
//         {
//             T = Ts[v][zi] ;
//             k = GET_ENCODING(v, zi) ;
//             Tdbg = ps[k] ;
//             if(T != Tdbg)
//             {
//                 ERR("T incorrect") ;
//             }
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Check locked arrays
//     /*------------------------------------------------------------------------*/
//     // Check that they are 0, 1
//     for(v = 0 ; v < 2 ; v++) Check_Array_Is_Binary(locked[v], n) ;

//     return ;
// }


// // Do all beneficial flips in either var 0, var 1, or both (depending on
// // whichVar). If withLocking is true, then lock the flips we do here.
// void Do_Beneficial_Flips
// (
//     VSPINT               *vars[2],
//     VSPINT              *Avars[2],
//     VSPINT               *wgts[2],
//     VSPINT             *locked[2],
//     VSPINT               Ts[2][2],
//     Doubly_LLs_Sized *Lists[2][2],
//     VSP_problem          *problem,
//     VSPINT              whichVars,
//     bool              withLocking
// )
// {
//     VSPINT indFlip, varFlip, gainFlip ;

//     for( ; 1 ; )
//     {
//         /*--------------------------------------------------------------------*/
//         // Select a flip
//         /*--------------------------------------------------------------------*/
//         Select_Flip
//         (
//             vars     , Avars    , wgts      , locked , Ts , Lists , problem,
//             &varFlip , &indFlip , &gainFlip, whichVars, FEAS_BEN_FLIP
//         ) ;

//         /*--------------------------------------------------------------------*/
//         // Do the best flip
//         /*--------------------------------------------------------------------*/
//         if(varFlip != NO_VAR)
//         {
//             Do_Flip
//             (
//                 vars    , Avars   , wgts    , locked, Ts, Lists,
//                 varFlip , indFlip , problem, withLocking
//             ) ;
//         }
//         else break ;
//     }
// }

// // Do upflips until the target weight is just about to be surpassed
// // ASSUMED:
// // - whichVar = 0 or 1, but not 2
// // - All downflips (so all zi = 1) are locked in whichVar
// //
// // NOT ASSUMED:
// // - Not assumed that all upflips are unlocked
// void Do_Upflips_Until_Target_Wgt
// (
//     VSPINT               *vars[2],
//     VSPINT              *Avars[2],
//     VSPINT               *wgts[2],
//     VSPINT             *locked[2],
//     VSPINT               Ts[2][2],
//     Doubly_LLs_Sized *Lists[2][2],
//     VSP_problem          *problem,
//     VSPINT               whichVar, // Can NOT be both variables
//     bool              withLocking,
//     Flip_Selector_Mode       mode,
//     VSPINT              targetWgt
// )
// {
//     VSPINT i, indFlip, gainFlip, varFlip,
//     n = problem->Graph->n,
//     *c = problem->Graph->c ;
//     /*------------------------------------------------------------------------*/
//     // Check inputs
//     /*------------------------------------------------------------------------*/
//     SET_FUNCTION_COUNTER ;
//     VSP_ASSERT(whichVar == 0 || whichVar == 1) ;

//     #ifdef VSP_DEBUG
//         // Check that all downflips are locked
//         for(i = 0 ; i < n ; i++)
//         {
//             if(vars[whichVar][i] && !locked[whichVar][i])
//             {
//                 ERR("Downflip not locked") ;
//                 break ;
//             }
//         }
//     #endif

//     // At the start of the loop, this is the number of flips done so far
//     VSPINT nflips = 0 ;
//     while( true )
//     {
//         // Select flip
//         Select_Flip
//         (
//             vars     , Avars    , wgts      , locked   , Ts , Lists , problem,
//             &varFlip , &indFlip , &gainFlip , whichVar , mode
//         ) ;

//         if(indFlip != EMPTY && *(wgts[whichVar]) + c[indFlip] <= targetWgt)
//         {
//             Do_Flip
//             (
//                 vars    , Avars   , wgts    , locked, Ts, Lists,
//                 varFlip , indFlip , problem, withLocking
//             ) ;
//             nflips++ ;
//         }
//         else break ;
//     }
//     return ;
// }
#include "shared_api.h"
// #include "mc_internal.h"

// #include <float.h>

// bool LP_is_optimal_case_1a(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) ;

// bool LP_is_optimal_case_1b(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) ;

// bool LP_is_optimal_case_2a(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) ;

// bool LP_is_optimal_case_2b(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) ;

// bool LP_is_optimal_case_3(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) ;

// /* Check that the given z is optimal for the LP
//     max (c - gamma Aw)' z
//     s.t. 0 <=  z   <= 1
//          1 <= c' z <= u
//     where u is the upper balance constraint
// */
// bool LP_is_optimal
// (
//     State         *Z,
//     Az           *Aw,
//     VSPFLOAT       gamma,
//     Problem *problem,
//     Options *options
// ) {

//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;

//         /* Check inputs */
//         check_state(Z, problem) ;

//         /* Compute the objective at the given Z, compute the optimal objective
//            by solving the LP, and compare results                             */

//         /* Compute objective at the given z                                   */
//         VSPFLOAT obj_given ;
//         obj_given = compute_LP_objective(Z->point, problem->Graph->c, gamma, Aw) ;

//         /* Get optimal objective by solving the LP                            */
//         VSPFLOAT obj_opt = 0 ;
//         State *Z_copy = Construct_State(problem->Graph->n) ;
//         copy_state(Z, Z_copy) ;
//         obj_opt = LP_solver_4(Z_copy      , Aw   , problem->Graph->c, 1,
//                                   problem->u, gamma, problem   , options, 'x');

//         UNUSED VSPFLOAT obj_opt2 = compute_LP_objective(Z_copy->point, problem->Graph->c, gamma, Aw) ;


//         /* Determine optimality                                               */
//         bool is_opt_dbg = equal(obj_given, obj_opt) ;

//         /* Free memory */
//         Destroy_State(Z_copy) ;
//     #endif

//     bool is_opt ;
//     /* Disjunct on the upper, lower, or neither balance contrains being active*/
//     if(Z->weight == problem->u) {
//         /* Disjunct on whether there is a fractional entry */
//         if(Z->point->ind_frac == EMPTY) {
//             is_opt = LP_is_optimal_case_1a(Z->point, Aw, gamma, problem->Graph->c);
//         }
//         else{
//             is_opt = LP_is_optimal_case_1b(Z->point, Aw, gamma, problem->Graph->c);
//         }
//     }
//     else if (Z->weight == 1) {
//         /* Disjunct on whether there is a fractional entry */
//         if(Z->point->ind_frac == EMPTY) {
//             is_opt = LP_is_optimal_case_2a(Z->point, Aw, gamma, problem->Graph->c);
//         }
//         else{
//             is_opt = LP_is_optimal_case_2b(Z->point, Aw, gamma, problem->Graph->c);
//         }
//     }
//     else {
//         is_opt = LP_is_optimal_case_3(Z->point, Aw, gamma, problem->Graph->c);
//     }

//     #ifdef DEBUG
//         VSPINT n = problem->Graph->n ;
//         (void) n ; // dummy
//         if(is_opt != is_opt_dbg)
//         {
//             Error("Error with optimality check", DBG_INFO) ;
//         }
//     #endif


//     return is_opt ;
// }




// /* c'z = u  and z is fully binary
//    Check that the maximum value density of the zeros is at most the minimum
//    value density of the 1's, and that the minimum value density of the ones
//    is non-negative.
//    In the code, this comes down to checking (Awi / ci) >= Awj / cj and that
//    cj - Awj >= 0 for i in the zeros and j in the ones                */
// bool LP_is_optimal_case_1a
// (
//     Point         *z,
//     Az           *Aw,
//     VSPFLOAT       gamma,
//     VSPINT            *c
// )
// {
//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;
//         if(z->ind_frac != EMPTY)
//         {
//             Error("Wrong case in optimality check", DBG_INFO) ;
//         }
//     #endif

//     VSPINT   n                 = Aw->n ;
//     VSPFLOAT Awi ;
//     VSPFLOAT val_frac_w        = Aw->val_frac_z ;

//     // unattainable lb on 1 - gamma Awi / ci is 1 - gamma * n - 1 = - gamma * n
//     // unattainable ub on 1 - gamma Awi / ci is 1
//     VSPFLOAT max_ratio_ones    = - gamma * n;
//     VSPFLOAT min_ratio_zeros   = 2 ; 
//     VSPFLOAT ratio             = 0 ;
//     VSPINT  *Aw_int            = Aw->int_part ;
//     VSPINT  *col_A             = Aw->col_frac_z ;
//     VSPINT  *z_int             = z->int_part ;
//     VSPINT   one               = Aw->one ;
//     VSPINT   w_has_ind_frac    = (Aw->ind_frac_z != EMPTY) ;

//     /* TODO speed: make different branches depending on the frac comp         */
//     for(VSPINT i = 0 ; i < n ; i++)
//     {    
//         Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ;
//         ratio = Awi / c[i] ;
//         if(z_int[i])
//         {
//             if(c[i] - gamma * Awi < 0) return false ;
//             else
//             {
//                 if(max_ratio_ones < ratio) max_ratio_ones = ratio ;
//             }
//         }
//         else
//         {
//             if(min_ratio_zeros > ratio) min_ratio_zeros = ratio ;
//         }

//         if (min_ratio_zeros < max_ratio_ones - ZERO) return false ;
//     }
//     return true ;
// }


// /* 1 < c'z = u and z has a fractional entry
//    Check that the max value density of the zeros bounded above by the value
//    density of the fractional component, which is itself a lower bound for the
//    min fractional density of the ones. Also, the value density of the
//    fractional entry must be non-negative. In the code, this means checking
//    Ayi / ci <= Ayk / ck <= Ayj / cj 
//    where i is an arbitrary index of the ones, k is the index of the fractional
//    entry, and j is an arbitrary index of the zeros, and that
//    ck - gamma * Ayk >= 0 */
// bool LP_is_optimal_case_1b(
//     Point         *z,
//     Az           *Aw,
//     VSPFLOAT       gamma,
//     VSPINT             *c
// ) {
//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;

//         if(z->ind_frac == EMPTY)
//         {
//             Error("Wrong case in optimality check", DBG_INFO) ;
//         }
//     #endif


//     VSPFLOAT Awi ;
//     VSPFLOAT val_frac_w        = Aw->val_frac_z ;
//     VSPINT  *Aw_int            = Aw->int_part ;
//     VSPINT  *col_A             = Aw->col_frac_z ;
//     VSPINT  *z_int             = z->int_part ;
//     VSPINT   one               = Aw->one ;
//     VSPINT   n                 = Aw->n ;
//     VSPINT   k                 = z->ind_frac ;
//     VSPINT   w_has_ind_frac    = (Aw->ind_frac_z != EMPTY) ;


//     VSPFLOAT Awk = Aw_int[k] + w_has_ind_frac * (col_A[k] == one) * val_frac_w ;

//     /* Value of fractional component should be non-negative */
//     if( c[k] - gamma * Awk < 0 ) return false ;

//     VSPFLOAT ratio_k = Awk / c[k] ;
//     for(VSPINT i = 0 ; i < n ; i++) {
//         Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ;
//         VSPFLOAT ratio_i = Awi / c[i] ;

//         /* For optimality, if zi = 0, ratio_i should be <= ratio_k, and
//            if zi = 1, ratio_i should be >= ratio_k */
//         VSPINT pm1 = (1 - 2 * z_int[i]) ;

//         // if ( pm1 * ratio_i < pm1 * ratio_k - ZERO )
//         if ( !lteq( pm1 * ratio_k, pm1 * ratio_i) )
//         {
//             return false ;
//         }
//     }
//     return true ;
// }




// /* 1 = c'z < u  and z is fully binary
//    Check that the maximum value density of the zeros is at most the minimum
//    value density of the 1's, and that the maximum value density of the zeros
//    is non-positive.
//    In the code, this comes down to checking (Awi / ci) >= Awj / cj and that
//    ci - Awi <= 0 for i in the zeros and j in the ones                */
// bool LP_is_optimal_case_2a(
//     Point         *z,
//     Az           *Aw,
//     VSPFLOAT       gamma,
//     VSPINT            *c
// ) {
//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;
//         if(z->ind_frac != EMPTY)
//         {
//             Error("Wrong case in optimality check", DBG_INFO) ;
//         }
//     #endif

//     VSPFLOAT Awi ;
//     VSPFLOAT val_frac_w        = Aw->val_frac_z ;
//     VSPFLOAT max_ratio_ones    = -DBL_MAX ;
//     VSPFLOAT min_ratio_zeros   =  DBL_MAX ;
//     VSPFLOAT ratio             = 0 ;
//     VSPINT  *Aw_int            = Aw->int_part ;
//     VSPINT  *col_A             = Aw->col_frac_z ;
//     VSPINT  *z_int             = z->int_part ;
//     VSPINT   one               = Aw->one ;
//     VSPINT   n                 = Aw->n ;
//     VSPINT   w_has_ind_frac    = (Aw->ind_frac_z != EMPTY) ;

//     /* TODO speed: make different branches depending on the frac comp         */
//     for(VSPINT i = 0 ; i < n ; i++)
//     {
//         Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ;
//         ratio = Awi / c[i] ;
//         if(z_int[i]) {
//             if(max_ratio_ones < ratio) max_ratio_ones = ratio ;
//         }
//         else {
//             if(c[i] - gamma * Awi > 0) return false ;
//             else {
//                 if(min_ratio_zeros > ratio) min_ratio_zeros = ratio ;
//             }
//         }
//         if (min_ratio_zeros < max_ratio_ones - ZERO) return false ;
//     }
//     return true ;
// }




// /* 1 = c'z < u and z has a fractional entry
//    Check that the max value density of the zeros bounded above by the value
//    density of the fractional component, which is itself a lower bound for the
//    min fractional density of the ones. Also, the value density of the
//    fractional entry must be non-positive. In the code, this means checking
//    Ayi / ci <= Ayk / ck <= Ayj / cj 
//    where i is an arbitrary index of the ones, k is the index of the fractional
//    entry, and j is an arbitrary index of the zeros, and that
//    ck - gamma Ayk <= 0 */
// bool LP_is_optimal_case_2b(
//     Point         *z,
//     Az           *Aw,
//     VSPFLOAT       gamma,
//     VSPINT             *c
// ) {
//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;
//         if(z->ind_frac == EMPTY)
//         {
//             Error("Wrong case in optimality check", DBG_INFO) ;
//         }
//     #endif

//     VSPFLOAT Awi ;
//     VSPFLOAT val_frac_w        = Aw->val_frac_z ;
//     VSPINT  *Aw_int            = Aw->int_part ;
//     VSPINT  *col_A             = Aw->col_frac_z ;
//     VSPINT  *z_int             = z->int_part ;
//     VSPINT   one               = Aw->one ;
//     VSPINT   n                 = Aw->n ;
//     VSPINT   k                 = z->ind_frac ;
//     VSPINT   w_has_ind_frac    = (Aw->ind_frac_z != EMPTY) ;

//     VSPFLOAT Awk = Aw_int[k] + w_has_ind_frac * (col_A[k] == one) * val_frac_w ;

//     /* Value of fractional component should be non-positive */
//     if( c[k] - gamma * Awk > 0 ) return false ;

//     VSPFLOAT ratio_k = Awk / c[k] ;
//     for(VSPINT i = 0 ; i < n ; i++) {
//         Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ;
//         VSPFLOAT ratio_i = Awi / c[i] ;

//         /* For optimality, if zi = 0, ratio_i should be <= ratio_k, and
//            if zi = 1, ratio_i should be >= ratio_k */
//         VSPINT pm1 = (1 - 2 * z_int[i]) ;
//         if ( pm1 * ratio_i < pm1 * ratio_k - ZERO )
//         {
//             return false ;
//         }
//     }
//     return true ;
// }





// /* z is fully binary and 1 < c' z < u.
//     Check:
//     for all i in Z: ci - gamma Awi <= 0
//     for all i in O: ci - gamma Awi >= 0  */
// bool LP_is_optimal_case_3(
//     Point          *z,
//     Az            *Aw,
//     VSPFLOAT        gamma,
//     VSPINT             *c
// ) {
//     #ifdef DEBUG
//         static long long unsigned count = 0 ;
//         count++ ;
//         if(z->ind_frac != EMPTY)
//         {
//             Error("Wrong case in optimality check", DBG_INFO) ;
//         }
//     #endif

//     VSPFLOAT Awi ;
//     VSPFLOAT val_frac_w     = Aw->val_frac_z ;
//     VSPINT  *Aw_int         = Aw->int_part ;
//     VSPINT  *col_A          = Aw->col_frac_z ;
//     VSPINT  *z_int          = z->int_part ;
//     VSPINT   one            = Aw->one ;
//     VSPINT   n              = Aw->n ;
//     VSPINT   w_has_ind_frac = (Aw->ind_frac_z != EMPTY) ;

//     for(VSPINT i = 0 ; i < n ; i++) {
//         /* TODO speed: perhaps splitting into cases is faster */
//         Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ; 
//         if(z_int[i]) {
//             /* TODO speed: profile and check if doing a < b instead is faster */
//             if(c[i] - gamma * Awi < -ZERO)
//             // if(c[i] < gamma * Awi )
//             {
//                 return false ;
//             }
//         }
//         else {
//             if (c[i] - gamma * Awi > ZERO) {
//                 return false ;
//             } 
//         }
//     }
//     return true ;
// }
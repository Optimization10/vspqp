#include "core.h"
#include "utils.h"
/*
   _________________________________________________________________________
  |                                                                         |
  | Check optimality of a feasible point of the Continuous Knapsack Problem |
  | arising in VSPLP.                                                       |
  |    max f(z) = v'z   subject to   c'z <= u ,  0 <= z <= 1, z in R^n      |
  | where v = c - gamma Aw                                                  |
  |                                                                         |
  | Assumed:                                                                |
  | 1) 0 < u < sum of all c_i                                               |
  | 2) n > 1                                                                |
  | 3) z is almost binary                                                   |
  | 4) c >= 1                                                               |
  |                                                                         |
  | For the theory, see the doc/Optimality_Conds_LP_No_lb.pdf               |
  |_________________________________________________________________________|

*/

typedef enum
{
    LP_OPT_CASE_1a1,
    LP_OPT_CASE_1a2,
    LP_OPT_CASE_1a3,
    LP_OPT_CASE_1a4,
    LP_OPT_CASE_1b1,
    LP_OPT_CASE_1b2,
    LP_OPT_CASE_1b3,
    LP_OPT_CASE_1b4,
    LP_OPT_CASE_2a1,
    LP_OPT_CASE_2a2,
    LP_OPT_CASE_2a3,
    LP_OPT_CASE_2a4,
    LP_OPT_CASE_2b1,
    LP_OPT_CASE_2b2,
    LP_OPT_CASE_2b3,
    LP_OPT_CASE_2b4
} LP_Opt_Checker_Case;


// Compute Value Density Ratio
static VSPFLOAT Val_Density(VSPINT i, VSPINT ci, VSPFLOAT gamma, Az *Aw)
{
    VSPFLOAT Awi = Azi(i, Aw) ;
    return 1 - Awi * gamma / ci ;
}


/* ========================================================================== */
// Cases 1a + 1,2,3,4
// - c'z = u
// - z has no fractional entry
/* ========================================================================== */

// Case 1a1:
// - c'z = u
// - z has no fractional entry
// - z has at least one 0 and at least one 1
// Checks:
// 1. min {r_j : z_j = 1} >= 0 (min val density of 1s is >= 0)
// 2. max {r_i : z_i = 0} <= min {r_j : z_j = 1} (max value density of 0s
//    is <= min value density of 1s)
static bool Is_LP_Optimal_Case_1a1
(
  VSPINT            *z,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;
    VSPFLOAT max_ri_zeros = - gamma * n ; // unattainable lower bnd on val dens
    VSPFLOAT min_ri_ones  = 2 ;           // unattainable upper bnd on val dens

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(z[i])
        {
            if( Is_Less_Than(ri, 0, tol) )
            {
                return false ;
            }
            if( ri < min_ri_ones )
            {
                min_ri_ones = ri ;
                if(Is_Greater_Than(max_ri_zeros, min_ri_ones, tol))
                {
                    return false ;
                }
            }
        }
        else // zi = 0
        {
            if( ri > max_ri_zeros)
            {
                max_ri_zeros = ri ;
                if(Is_Greater_Than(max_ri_zeros, min_ri_ones, tol))
                {
                    return false ;
                }
            }
        }
    }
    return true ;
}


// Case 1a2:
// - c'z = u
// - z has no fractional entry
// - z has at least one 0 and no 1s (so z is all 0s since no frac entry)
//
// This case is impossible because z is all 0s -> c'z = 0 = u, but u is assumed
// to be > 0

// Case 1a3:
// - c'z = u
// - z has no fractional entry
// - z has no 0s but it has at least one 1 (so z is all 1s since no frac entry)
//
// This case is impossible because z is all 1s -> z'c = 1'c, and this 
// together with c'z = u implies u = 1'c which contradicts our assumptions


// Case 1a4:
// - c'z = u
// - z has no fractional entry
// - z has no 0s and no 1s
//
// This case is impossible since z has no fractional entry either

/* ========================================================================== */
/* ========================================================================== */
// Cases 1b + 1,2,3,4
// - c'z = u
// - z has a fractional entry in index k  
/* ========================================================================== */
/* ========================================================================== */

// Case 1b1:
// - c'z = u
// - z has a fractional entry in index k
// - z has at least one 0 and at least one 1
//
// Checks:
// 1. rk >= 0
// 2. max {r_i : z_i = 0} <= r_k
// 3. r_k <= min {r_j : z_j = 1}
static bool Is_LP_Optimal_Case_1b1
(
    VSPINT            *z,
    VSPINT             k,
    VSPFLOAT          zk,
    Az           *Aw,
    VSPINT            *c,
    VSPINT             u,
    VSPFLOAT       gamma,
    VSPFLOAT         tol
)
{
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if(Is_Less_Than(rk, 0, tol)) return  false ;

    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(z[i])
        {
            if(Is_Less_Than(ri, rk, tol))
            {
                return false ;
            }
        }
        else
        {
            if(Is_Greater_Than(ri, rk, tol))
            {
                return false ;
            }
        }
    }
    return true ;
}

// Case 1b1:
// - c'z = u
// - z has a fractional entry in index k
// - z has at least one 0 and no 1s
// Checks:
// 1. rk >= 0
// 2. max {r_i : z_i = 0} <= r_k
static bool Is_LP_Optimal_Case_1b2
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if(Is_Less_Than(rk, 0, tol)) return  false ;

    VSPINT i, n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        #ifdef DEBUG
            VSP_ASSERT(z[i] == 0) ;
        #endif

        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(Is_Greater_Than(ri, rk, tol))
        {
            return false ;
        }
    }
    return true ;
}


// Case 1b3:
// - c'z = u
// - z has a fractional entry in index k
// - z has no 0s and it has at least one 1
// Checks:
// 1. rk >= 0
// 2. min {r_j : z_j = 1} => r_k
static bool Is_LP_Optimal_Case_1b3
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if(Is_Less_Than(rk, 0, tol)) return  false ;

    VSPINT i, n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        #ifdef DEBUG
            if(z[i] != 1 && i != k)
            {
                ERR("Wrong opt checker case") ;
                break ;
            }
        #endif

        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(Is_Less_Than(ri, rk, tol))
        {
            return false ;
        }
    }
    return true ;
}


// Case 1b4:
// - c'z = u
// - z has a fractional entry in index k
// - z has no 0s and no 1s
//
// This case is not possible because n > 1 implies there is at least one 0
// or one 1



/* ========================================================================== */
/* ========================================================================== */
// Cases 2a + 1,2,3,4
// - c'z < u
// - z has no fractional entry  
/* ========================================================================== */
/* ========================================================================== */
// Case 2a1:
// - c'z < u
// - z has no fractional entry
// - z has at least one 0 and at least one 1
//
// Check:
// 1. r_i <=  0 for all i with z_i = 0
// 1. r_i >=  0 for all i with z_i = 1
static bool Is_LP_Optimal_Case_2a1
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(z[i])
        {
            if(Is_Less_Than(ri, 0, tol))
            {
                return false ;
            }
        }
        else
        {
            if(Is_Greater_Than(ri, 0, tol))
            {
                return false ;
            }
        }
    }
    return true ;
}


// Case 2a2:
// - c'z < u
// - z has no fractional entry
// - z has at least one 0 and no 1s
//
// Check:
// 1. r_i <=  0 for all i with z_i = 0
static bool Is_LP_Optimal_Case_2a2
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        VSP_ASSERT(z[i] == 0) ;
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(Is_Greater_Than(ri, 0, tol))
        {
            return false ;
        }
    }
    return true ;
}


// Case 2a3:
// - c'z < u
// - z has no fractional entry
// - z has no 0s and at least one 1
//
// This case is not possible since z is all 1s, so c'z = c'1 > u


// Case 2a4:
// - c'z < u
// - z has no fractional entry
// - z has no 0s and no 1s
//
// This case is not possible since n > 1 and z is almost binary


/* ========================================================================== */
/* ========================================================================== */
// Cases 2b + 1,2,3,4
// - c'z < u
// - z has a fractional entry in index k  
/* ========================================================================== */
/* ========================================================================== */
// Case 2b1:
// - c'z < u
// - z has a fractional entry in index k
// - z has at least one 0 and at least one 1
//
// Check:
// 1. rk = 0
// 2. max {r_i : z_i = 0} <= 0
// 3. 0 <= min {r_j : z_j = 1}

static bool Is_LP_Optimal_Case_2b1
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if( !equal(rk, 0) ) return  false ;

    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(z[i])
        {
            if(Is_Less_Than(ri, 0, tol))
            {
                return false ;
            }
        }
        else
        {
            if(Is_Greater_Than(ri, 0, tol))
            {
                return false ;
            }
        }
    }
    return true ;
}


// Case 2b1:
// - c'z < u
// - z has a fractional entry in index k
// - z has at least one 0 and no 1s
//
// Check:
// 1. rk = 0
// 2. max {r_i : z_i = 0} <= 0
static bool Is_LP_Optimal_Case_2b2
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if( !equal(rk, 0) ) return  false ;

    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        VSP_ASSERT(z[i] == 0) ;
        if(Is_Greater_Than(ri, 0, tol))
        {
            return false ;
        }
    }
    return true ;
}

// Case 2b3:
// - c'z < u
// - z has a fractional entry in index k
// - z has no 0s and at least one 1
//
// Check:
// 1. rk = 0
// 2. r_i >= 0 for all i
static bool Is_LP_Optimal_Case_2b3
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;
    VSPFLOAT rk = Val_Density(k, c[k], gamma, Aw) ;
    if( !equal(rk, 0) ) return  false ;

    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;

    #ifdef DEBUG
        for(i = 0 ; i < n ; i++)
        {
            if(z[i] != 1 && i != k)
            {
                ERR("Wrong opt checker case") ;
                break ;
            }
        }
    #endif

    for(i = 0 ; i < n ; i++)
    {
        ri = Val_Density(i, c[i], gamma, Aw) ;
        VSP_ASSERT(z[i] == 1) ;
        if(Is_Less_Than(ri, 0, tol))
        {
            return false ;
        }
    }
    return true ;
}


// Case 2b4:
// - c'z < u
// - z has a fractional entry in index k
// - z has no 0s and no 1s
//
// Not possible since z is almost binary and n > 1

/* ========================================================================== */
/* ========================================================================== */
// Optimality Checker Manager
/* ========================================================================== */
/* ========================================================================== */


static void check_has0s_and_has1s
(
    VSPINT        *z,
    VSPINT         n,
    VSPINT         k,
    bool has0s_claim,
    bool has1s_claim
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    #ifndef DEBUG
        return ;
    #endif

    bool has0s = false ;
    bool has1s = false ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        if(i == k) continue ; // skip frac entry
        if(z[i]) has1s = true ;
        else     has0s = true ;
    }

    if(has0s != has0s_claim)
    {
        ERR("LP Opt claim mismatch") ;
    }

    if(has1s != has1s_claim)
    {
        ERR("LP Opt claim mismatch") ;
    }
    return ;
}


// z has a 1 in each entry
bool All_Ones_Case
(
  VSPINT            *z,
  VSPINT             k,
  VSPFLOAT          zk,
  Az           *Aw,
  VSPINT            *c,
  VSPINT             u,
  VSPFLOAT       gamma,
  VSPFLOAT         tol
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(k == EMPTY) ;

    VSPINT i, j,
    n = Aw->n ;
    VSPFLOAT ri ;

    for(i = 0 ; i < n ; i++)
    {
        VSP_ASSERT(z[i] == 1) ;
        ri = Val_Density(i, c[i], gamma, Aw) ;
        if(Is_Less_Than(ri, 0, tol))
        {
            return false ;
        }
    }
    return true ;
}

bool Is_LP_Optimal
(
  VSPINT              *z,
  VSPINT               k,
  VSPFLOAT            zk,
  VSPINT           wgt_z,
  Az                 *Aw,
  VSPINT              *c,
  VSPINT               u,
  VSPFLOAT         gamma,
  VSPINT      sum_all_ci,
  Options       *options
)
{
    SET_FUNCTION_COUNTER ;
    bool ub_active      = (wgt_z == u) ;
    bool has_frac_entry = (k != EMPTY) ;
    bool has_0s ;
    bool has_1s ;
    VSPFLOAT tol = options->tol_LP_opt_checker ;

    if(has_frac_entry)
    {
        VSPFLOAT wgt_frac            = zk * c[k] ;
        has_1s = ( Is_Greater_Than(wgt_z, wgt_frac, tol) ) ;

        // if z is all 1s, then rounding up the fractional entry would result
        // in c'z = sum of all ci. Conversely, if c'z is not summ of all c_i,
        // then z has at least one 0
        VSPFLOAT wgt_frac_complement = (1 - zk) * c[k] ;
        has_0s = Is_Less_Than(wgt_z + wgt_frac_complement, sum_all_ci, tol) ;
    }
    else // no frac entry
    {
        // Should not happen
        has_0s = (wgt_z < sum_all_ci) ;
        has_1s = (wgt_z > 0) ;
    }

    #ifdef DEBUG
        check_has0s_and_has1s(z, Aw->n, k, has_0s, has_1s) ;
    #endif

    LP_Opt_Checker_Case Case ;
    bool case_1       = ub_active          ;
    bool subcase_a    = !has_frac_entry    ;
    bool subsubcase_1 = has_0s  && has_1s  ;
    bool subsubcase_2 = has_0s  && !has_1s ;
    bool subsubcase_3 = !has_0s && has_1s  ;
    bool subsubcase_4 = !has_0s && !has_1s ;

    if(case_1)
    {
        if(subcase_a)
        {
            if     (subsubcase_1) Case = LP_OPT_CASE_1a1 ;
            else if(subsubcase_2) Case = LP_OPT_CASE_1a2 ;
            else if(subsubcase_3) Case = LP_OPT_CASE_1a3 ;
            else                  Case = LP_OPT_CASE_1a4 ;
        }
        else // subcase_b
        {
            if     (subsubcase_1) Case = LP_OPT_CASE_1b1 ;
            else if(subsubcase_2) Case = LP_OPT_CASE_1b2 ;
            else if(subsubcase_3) Case = LP_OPT_CASE_1b3 ;
            else                  Case = LP_OPT_CASE_1b4 ;
        }
    }
    else // case_2
    {
        if(subcase_a)
        {
            if     (subsubcase_1) Case = LP_OPT_CASE_2a1 ;
            else if(subsubcase_2) Case = LP_OPT_CASE_2a2 ;
            else if(subsubcase_3) Case = LP_OPT_CASE_2a3 ;
            else                  Case = LP_OPT_CASE_2a4 ;
        }
        else // subcase_b
        {
            if     (subsubcase_1) Case = LP_OPT_CASE_2b1 ;
            else if(subsubcase_2) Case = LP_OPT_CASE_2b2 ;
            else if(subsubcase_3) Case = LP_OPT_CASE_2b3 ;
            else                  Case = LP_OPT_CASE_2b4 ;
        }
    }

    // Bandaid: Short-circuir the case of all 1s, which is temporarily possible
    // because I haven't removed the compression yet
    if(k == EMPTY && !has_0s)
    {
        return All_Ones_Case(z, k, zk, Aw, c, u, gamma, tol) ;
    }

    switch (Case)
    {
        case LP_OPT_CASE_1a1:
            return Is_LP_Optimal_Case_1a1(z, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_1a2:
            printf("\nError: Impossible case 1a2 reached\n");
            return false;
        case LP_OPT_CASE_1a3:
            printf("\nError: Impossible case 1a3 reached\n");
            return false;
        case LP_OPT_CASE_1a4:
            printf("\nError: Impossible case 1a4 reached\n");
            return false;
        case LP_OPT_CASE_1b1:
            return Is_LP_Optimal_Case_1b1(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_1b2:
            return Is_LP_Optimal_Case_1b2(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_1b3:
            return Is_LP_Optimal_Case_1b3(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_1b4:
            printf("\nError: Impossible case 1b4 reached\n");
            return false;
        case LP_OPT_CASE_2a1:
            return Is_LP_Optimal_Case_2a1(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_2a2:
            return Is_LP_Optimal_Case_2a2(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_2a3:
            printf("\nError: Impossible case 2a3 reached\n");
            return false;
        case LP_OPT_CASE_2a4:
            printf("\nError: Impossible case 2a4 reached\n");
            return false;
        case LP_OPT_CASE_2b1:
            return Is_LP_Optimal_Case_2b1(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_2b2:
            return Is_LP_Optimal_Case_2b2(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_2b3:
            return Is_LP_Optimal_Case_2b3(z, k, zk, Aw, c, u, gamma, tol);
        case LP_OPT_CASE_2b4:
            printf("\nError: Impossible case 2b4 reached\n");
            return false;
        default:
            printf("\nError: Case not handled by LP opt checker") ;
            return false ;
    }

    printf("\nError: Reached impossible point") ;
    return false; // unreachable
}
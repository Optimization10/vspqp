#include "shared_api.h"
#include "mc_internal.h"
#include "ds_api.h"
#include <float.h>

/*
* Solve knapsack problem
*         max  v * z
*         s.t.       c * z  <= u,
*               0 <=   z    <= 1
* where v = c - gamma Aw, efficiently.
*
* The algorithm is as follows:
* Partition 1, 2, 3, ..., n into
*     L1 = {i : v_i / c_i = 1       } = {i : Aw_i = 0               }
*     L2 = {i : v_i / c_i in (0, 1) } = {i : 0 < gamma * Aw_i < c_i }
*     L3 = {i : v_i / c_i <= 0      } = {i : ci <= gamma * Aw_i     }
*/

/* Does NOT use partial max sort */
VSPFLOAT LP_solver_4
(
    State              *Z, /* Solution */
    Az                *Aw,
    VSPINT                 *c,
    VSPINT                  u,
    VSPFLOAT            gamma,
    Problem      *problem,
    Options      *options,
    char                  var  /* Variable we are optimizing over */
)
{

    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        VSPFLOAT objStart = compute_LP_objective(Z->point, c, gamma, Aw) ;
    #endif


    Az *Az = Z->Az ;
    VSPINT n = problem->Graph->n ;
    VSPINT z_weight = 0 ;
    VSPINT i, ci, k ;
    VSPFLOAT Awi ;
    VSPINT *Aw_int = Aw->int_part ;
    VSPINT const ind_frac_w = Aw->ind_frac_z ;
    VSPFLOAT val_frac_w = Aw->val_frac_z ;
    VSPINT *col_A = Aw->col_frac_z ;
    VSPINT one = Aw->one ;

    VSPINT *z_int = Z->point->int_part ;
    VSPINT *L2 = malloc( n * sizeof(*L2) ) ;
    VSPINT n_L2 = 0 ;
    VSPINT L2_weight = 0 ;
    VSPFLOAT obj = 0 ;

    // VSPINT min_L2_weight = problem->Graph->max_vtx_cost ;

    /* Reset z fractional index */
    Z->point->ind_frac = EMPTY ;
    Z->Az->ind_frac_z  = EMPTY ;

    bool useRandomness = false ;
    // bool useRandomness = (Z->weight == problem->u) ;
    VSPINT *perm ;
    VSPINT j ;
    if(useRandomness)
    {
        perm = malloc(n * sizeof(VSPINT)) ;
        Get_Permutation(perm, n) ;
    }

    bool reverse = true ;
    if(rand()%2)
    {
        reverse = true ;
    }

    for(j = 0 ; j < n ; j++)
    {
        
        i = (useRandomness) ? perm[j] : j ;

        if(reverse) i = n - 1 - j ;
        else i = j ;

        #ifdef DEBUG
            VSPINT weight_dbg = 0 ;
            for(VSPINT p = 0 ; p < i ; p++)
            {
                weight_dbg += z_int[p] * c[p] ;
            }
        #endif


        /* TODO: Try branchless code here and see if it makes a difference */
        if(ind_frac_w == EMPTY)
        {
            Awi = Aw_int[i] ;
        }
        else
        {
            Awi = Aw_int[i] + (col_A[i] == one) * val_frac_w ;
        }
        ci = c[i] ;

        /* i is in L1 */
        if(Awi == 0) 
        {
            /* Set zi to 1 and continue */
            if(z_weight + ci < u)
            {
                z_weight += ci ;
                obj      += ci - gamma * Awi ;
                if(z_int[i] != 1)
                {
                    z_int[i] = 1 ;
                    Az_pm_col_of_A(Az->int_part, i, PLUS, problem) ;
                }
                continue ;
            }
            /* Set zi to 1, finish up and return */
            else if(z_weight + ci == u)
            {
                obj      += ci - gamma * Awi ;
                if(z_int[i] != 1)
                {
                    z_int[i] = 1 ;
                    Az_pm_col_of_A(Az->int_part, i, PLUS, problem) ;
                }

                /* Finish up by setting all remaining i to 0 */
                for(j++ ; j < n ; j++)
                {
                    i = (useRandomness) ? perm[j] : j ;
                    if(reverse) i = n - 1 - j ;
                    else i = j ;

                    if(z_int[i] != 0)
                    {
                        z_int[i] = 0 ;
                        Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
                    }
                }
                Z->weight = u ;
                free(L2) ;

                #ifdef DEBUG
                    check_state(Z, problem) ;
                    /* Check KKT conditions */
                    bool is_opt = Is_LP_Optimal
                    (
                        Z->point->int_part,
                        Z->point->ind_frac,
                        Z->point->val_frac,
                        Z->weight,
                        Aw,
                        c, u, gamma, problem->Graph->total_vtx_cost,
                        options
                    ) ;

                    if(is_opt == false) 
                    {
                        Error("Issue in LP solver", DBG_INFO) ;
                    }
                    VSPFLOAT objEnd = compute_LP_objective(Z->point, c, gamma, Aw) ;
                    VSP_ASSERT(objStart <= objEnd) ;
                #endif
                if(useRandomness) free(perm) ;
                return obj ;
            }
            /* If z_weight +ci > u. Set zi to a fraction, finish up and return*/
            else
            {
                Z->point->ind_frac = i ;
                Z->point->val_frac = (u - z_weight) / ((VSPFLOAT) ci) ;
                Z->Az->ind_frac_z  = i ;
                Z->Az->val_frac_z  = (u - z_weight) / ((VSPFLOAT) ci) ;
                obj               += (ci - gamma * Awi) * Z->point->val_frac ;

                Z->Az->one++ ;
                extract_matrix_col
                (
                    Z->Az->col_frac_z,
                    Z->point->ind_frac,
                    Z->Az->one,
                    problem->Graph->Ai,
                    problem->Graph->Ap
                ) ;

                /* Finish up by setting all remaining i to 0. Start with
                   current i, because it is frac index */
                for( ; j < n ; j++)
                {
                    i = (useRandomness) ? perm[j] : j ;
                    if(reverse) i = n - 1 - j ;
                    else i = j ;

                    if(z_int[i] != 0)
                    {
                        z_int[i] = 0 ;
                        Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
                    }
                }
                Z->weight = u ;
                free(L2) ;
                #ifdef DEBUG
                    check_state(Z, problem) ;

                    /* Check KKT conditions */
                    bool is_opt = Is_LP_Optimal
                    (
                        Z->point->int_part,
                        Z->point->ind_frac,
                        Z->point->val_frac,
                        Z->weight,
                        Aw,
                        c, u, gamma, problem->Graph->total_vtx_cost,
                        options
                    ) ;

                    if(is_opt == false) 
                    {
                        Error("Issue in LP solver", DBG_INFO) ;
                    }

                    VSPFLOAT objEnd = compute_LP_objective(Z->point, c, gamma, Aw) ;
                    VSP_ASSERT(objStart <= objEnd) ;
                #endif
                if(useRandomness) free(perm) ;
                return obj ;
            }
        }
        /* i is in L3 */
        else if(ci <= Awi * gamma) /* TODO misc: consider +/- 1E-10 here */
        // else if(ci <= Awi * gamma - 1E-8) /* This seems to give better results, but it puts indices into L2 */
        {
            /* Set zi to 0 */
            if(z_int[i] != 0)
            {
                z_int[i] = 0 ;
                Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
            }
        }
        /* i is in L2 */
        else
        {
            /* Add i to L2. Compute the weight of L2. Set zi to 0  */
            /* TODO speed: consider if to set these to 0 or not. If not,
               then make sure than in the early exit above, we set them to 0
               before returning */
            L2[ n_L2++ ] = i ;
            L2_weight += ci ;

            // if(ci < min_L2_weight)
            // {
            //     min_L2_weight = ci ;
            // }
            // TODO: I think there should be a better way to do this
            if(z_int[i] != 0)
            {
                z_int[i] = 0 ;
                Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
            }
        }
    }

    /* If all of L2 fits into the knapsack, take it and return */
    if(z_weight + L2_weight <= u)
    {
        /* Take all items in L2 and return */
        Z->weight = z_weight + L2_weight ;
        for(k = 0 ; k < n_L2 ; k++)
        {
            i = L2[k] ;
            if(ind_frac_w == EMPTY)
            {
                Awi = Aw_int[i] ;
            }
            else
            {
                Awi = Aw_int[i] + (col_A[i] == one) * val_frac_w ;
            }

            ci = c[i] ;
            obj += ci - gamma * Awi ;
            if(z_int[i] != 1)
            {
                z_int[i] = 1 ;
                Az_pm_col_of_A(Az->int_part, i, PLUS, problem) ;
            }
        }
        free(L2) ;
        #ifdef DEBUG
            check_state(Z, problem) ;

            /* Check KKT conditions */
            bool is_opt = Is_LP_Optimal
            (
                Z->point->int_part,
                Z->point->ind_frac,
                Z->point->val_frac,
                Z->weight,
                Aw,
                c, u, gamma, problem->Graph->total_vtx_cost,
                options
            ) ;
            if(is_opt == false) 
            {
                Error("Issue in LP solver", DBG_INFO) ;
            }
            VSPFLOAT objEnd = compute_LP_objective(Z->point, c, gamma, Aw) ;
            VSP_ASSERT(objStart <= objEnd) ;
        #endif

        if(useRandomness) free(perm) ;
        return obj ;
    }

    /* In this case, we repeatedly take the minimal element of L2 until the
       knapsack is full */
    min_heap_build_v3(L2 - 1, Aw, c, n_L2 ) ;

    while(true)
    {
        VSP_ASSERT(n_L2 > 0) ;
        i = L2[0] ;
        ci = c[i] ;
        if(ind_frac_w == EMPTY)
        {
            Awi = Aw_int[i] ;
        }
        else
        {
            Awi = Aw_int[i] + (col_A[i] == one) * val_frac_w ;
        }

        if(z_weight + ci < u)
        {
            /* Take i fully */
            z_weight += ci ;
            obj += ci - gamma * Awi ;
            if(z_int[i] != 1)
            {
                z_int[i] = 1 ;
                Az_pm_col_of_A(Az->int_part, i, PLUS, problem) ;
            }
            min_heap_delete_v3(L2 - 1, Aw, c, &n_L2) ;
            continue ;
        }
        else if(z_weight + ci == u)
        {
            /* Take i fully */
            z_weight += ci ;
            obj += ci - gamma * Awi ;
            if(z_int[i] != 1)
            {
                z_int[i] = 1 ;
                Az_pm_col_of_A(Az->int_part, i, PLUS, problem) ;
            }
            break ;
        }
        else
        {
            /* Take a fractional part of i */
            Z->point->ind_frac = i ;
            Z->point->val_frac = (u - z_weight) / ((VSPFLOAT) ci) ;
            Z->Az->ind_frac_z  = i ;
            Z->Az->val_frac_z  = (u - z_weight) / ((VSPFLOAT) ci) ;
            obj               += (ci - gamma * Awi) * Z->point->val_frac ;

            if(z_int[i] != 0)
            {
                z_int[i] = 0 ;
                Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
            }

            Z->Az->one++ ;
            extract_matrix_col
            (
                Z->Az->col_frac_z,
                Z->point->ind_frac,
                Z->Az->one,
                problem->Graph->Ai,
                problem->Graph->Ap
            ) ;
            break ;
        }
    }

    /* Set any remaining members of L2 to 0 */
    for(k = 1 ; k < n_L2 ; k++)
    {
        i = L2[k] ;
        if(z_int[i] != 0)
        {
            z_int[i] = 0 ;
            Az_pm_col_of_A(Az->int_part, i, MINUS, problem) ;
        }
    }

    Z->weight = u ;
    free(L2) ;

    #ifdef DEBUG
        check_state(Z, problem) ;

        /* Check KKT conditions */
        bool is_opt = Is_LP_Optimal
        (
            Z->point->int_part,
            Z->point->ind_frac,
            Z->point->val_frac,
            Z->weight,
            Aw,
            c, u, gamma, problem->Graph->total_vtx_cost,
            options
        ) ;

        if(is_opt == false) 
        {
            Error("Issue in LP solver", DBG_INFO) ;
        }
        VSPFLOAT objEnd = compute_LP_objective(Z->point, c, gamma, Aw) ;
        VSP_ASSERT(objStart - 1E-8 <= objEnd) ;
        if( !(objStart - 1E-8 <= objEnd ) )
        {
            printf("\n\nobj start: %f \nobj end : %f", objStart, objEnd) ;
        }
    #endif
    return obj ;

}




/* ============================================================================
   ================ LP_solver_0 ===========================================
   ============================================================================

 PURPOSE: Solve knapsack problem
        max  d * x
        s.t. l  <= c * x <= u
             0  <= x     <= 1
 by sorting all the value densities
*/

VSPFLOAT LP_solver_0
(
    VSPFLOAT  *x,         /* solution                                         */
    VSPFLOAT  *d,         /* objective                                        */
    VSPINT    *c,           /* constraint                                     */
    VSPINT     l,            /* Lower bound. Must be 0 or 1                   */
    VSPINT     u,            /* RHS of constraint. Must be at least 1         */
    VSPINT     n,            /* size of problem                               */
    Options *options /* Choose randomness or not                          */
)
{
#ifdef DEBUG
    static long long unsigned count = 0;
    count += 1;
#endif

    VSPINT i = 0, weight = 0;
    VSPFLOAT obj_val;
    /* Calculate the value density array r */
    VSPFLOAT *r = calloc(n, sizeof(VSPFLOAT));
    for (i = 0; i < n; i++) r[i] = d[i] / c[i];

    /* store r in a bunch of structs, so we can sort it */
    value_index *r_structs = calloc(n, sizeof(value_index));
    for (i = 0; i < n; i++)
    {
        r_structs[i].float_value = r[i];
        r_structs[i].index = i;
    }

    /* sort the structs */
    qsort(r_structs, n, sizeof(value_index), compare_value_index_floats);

    /* solve the problem */
    for (i = 0; i < n; i++) x[i] = 0;
    VSPINT take_fractional_item = 0;
    VSPINT index = 0;

    /* Case if lower bound is 1 */
    if (l == 1 && r_structs[0].float_value <= 0)
    {
        /* Take as little as possible of the least worst item, and return */
        i = r_structs[0].index;
        x[i] = 1 / (double)c[i];

        // calculate objective value
        obj_val = 0;
        for (i = 0; i < n; i++)
        {
            obj_val += x[i] * d[i];
        }

        /* Free memory and return */
        free(r);
        free(r_structs);
        return obj_val;
    }

    /* Case when lower bound is 0 */
    for (i = 0; i < n; i++)
    {
        value_index entry = r_structs[i];
        index = entry.index;
        if (entry.float_value > 0)
        {
            if (weight + c[index] <= u)
            {
                x[index] = 1;
                weight += c[index];
            }
            else
            {
                take_fractional_item = 1;
                break;
            }
        }
    }
    if (take_fractional_item)
    {
        x[index] = (VSPFLOAT)(u - weight) / c[index];
        weight = u;
    }

    /* Calculate objective value */
    obj_val = 0;
    for (i = 0; i < n; i++)
    {
        obj_val += x[i] * d[i];
    }

    #ifdef DEBUG
        // test if r is ordered in decreasing order
        for (i = 0; i < n - 1; i++)
        {
            if ((r_structs[i]).float_value < (r_structs[i + 1]).float_value)
            {
                Error("Value densities not sorted correctly", DBG_INFO);
                break ;
            }
        }


        // make sure all entries of x are 0 or 1
        for (i = 0; i < n; i++)
        {
            if (x[i] < 0 || x[i] > 1)
            {
                Error("Entry not binary", DBG_INFO);
                break ;
            }
        }

        // check if the solution is feasible
        VSPINT weight_check = 0;
        for (i = 0; i < n; i++) weight_check += x[i] * c[i];
        if (weight_check > u) Error("Solution not feasible", DBG_INFO);
    #endif

    // free memory
    free(r);
    free(r_structs);

    return obj_val;
}
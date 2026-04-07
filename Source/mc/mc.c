#include "core.h"
#include "shared_api.h"
#include "mc_internal.h"
#include "shared_functions.h"
#include "utils.h"
#include <stdio.h>

#define EPS (1E-10)

void mountain_climbing
(
    State               *X,
    State               *Y,
    Problem       *problem,
    Options       *options,
    VSPFLOAT             gamma,
    bool          is_gamma_org
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        if(is_gamma_org) VSP_ASSERT(equal(gamma,problem->graph->max_vtx_cost)) ;
        else             VSP_ASSERT(gamma < problem->graph->max_vtx_cost) ;
    #endif
    TIMER_START ;

    #ifdef DEBUG
        /* Check input */
        check_state(Y, problem);
        check_state(X, problem);

        VSPFLOAT obj_dbg ;
        obj_dbg = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    /* Local vars */
    VSPINT *c  = problem->graph->c;
    VSPINT bal = problem->u;
    VSPINT n   = problem->graph->n;

    State *Xstar = Clone_State(X) ;
    State *Ystar = Clone_State(Y);

    VSPFLOAT obj = compute_obj(X, Y, problem, gamma);
    VSPINT loop = 0;
    char optimal_var;

    /*** Diagonal Phase *******************************************************/
    while (true)
    {
        loop++;

        /* Optimize over x to get x_star */
        VSPFLOAT f_xstar_y = LP_solver_4(Xstar, Y->Az, c, bal,
                                                 gamma, problem, options, 'x') ;
        f_xstar_y += Y->weight;

        /* Optimize over y to get y_star */
        VSPFLOAT f_x_ystar = LP_solver_4(Ystar, X->Az, c, bal,
                                                 gamma, problem, options, 'y') ;
        f_x_ystar += X->weight;

        /* Adjust objective value to prevent cycling */
        f_xstar_y = MAX(f_xstar_y, obj);
        f_x_ystar = MAX(f_x_ystar, obj);

        /* Compare f(x*, y*) vs f(x_star, y) vs f(x, y_star) */
        VSPFLOAT f_xstar_ystar;
        f_xstar_ystar = compute_obj(Xstar, Ystar, problem, gamma);

        if (f_xstar_ystar > MAX(f_xstar_y, f_x_ystar) + MAX(f_xstar_y, f_x_ystar) * EPS)
        {
            /* Set (x, y) <- (x_star, y_star) and repeat */
            /* TODO speed: Don't have 2 states and copy. Just maintain 4 states
               throughout an entire level                                     */
            copy_state(Xstar, X);
            copy_state(Ystar, Y);

            obj = f_xstar_ystar;
            continue;
        }
        else if (f_xstar_y > f_x_ystar + f_x_ystar * EPS)
        {
            copy_state(Xstar, X);
            optimal_var = 'x';
            break ;
        }
        else
        {
            copy_state(Ystar, Y);
            optimal_var = 'y';
            break;
        }
    }

    /* Free memory */
    Destroy_State(Xstar);
    Destroy_State(Ystar);

    /*** Orthogonal Phase *****************************************************/
    VSPFLOAT old_obj = compute_obj(X, Y, problem, gamma);
    VSPINT mc_count = 0 ;
    while (true)
    {
        VSPFLOAT new_obj;

        if (optimal_var == 'x')
        {
            /* Optimize over y */
            new_obj = LP_solver_4(Y, X->Az, c, bal, gamma, problem, options, 'y') ;
            new_obj += X->weight;
            optimal_var = 'y';
        }
        else
        {
            /* Optimize over x */
            new_obj = LP_solver_4(X, Y->Az, c, bal, gamma, problem, options, 'x') ;
            new_obj += Y->weight;
            optimal_var = 'x';
        }

        mc_count++ ;
        if(mc_count > 3) break ;

        new_obj = MAX(old_obj, new_obj);

        /* TODO: Consider putting an option to disable this optimality checker.
           The cost in CPU time is relatively high for the benefit we get */
        if (new_obj <= old_obj + 1E-8)
        {
            bool is_opt ;
            if(optimal_var == 'y') // Since y was just optimized over, we care
                                   // if x is optimal
            {
                is_opt = Is_LP_Optimal
                (
                    X->point->int_part,
                    X->point->ind_frac,
                    X->point->val_frac,
                    X->weight,
                    Y->Az,
                    c, problem->u, gamma, problem->graph->total_vtx_cost,
                    options
                ) ;
                if(is_opt) break ;
            }
            else
            {
                is_opt = Is_LP_Optimal
                (
                    Y->point->int_part,
                    Y->point->ind_frac,
                    Y->point->val_frac,
                    Y->weight,
                    X->Az,
                    c, problem->u, gamma, problem->graph->total_vtx_cost,
                    options
                ) ;
                if(is_opt) break ;
            }
        }

        old_obj = new_obj;
    }
    TIMER_END ;

    // printf("level: %lld  n: %lld  nz: %lld gamma_org: %lld gamma_used: %f mc_count: %lld  mc_time: %Lf\n",
    // PRI(problem->level), PRI(n), PRI(problem->graph->nz), PRI(problem->graph->max_vtx_cost),
    // gamma, PRI(mc_count),
    // TIMER_ELAPSED ) ;

    // if(!is_gamma_org)
    // {
    //     printf("%s,%lld,%lld,%lld,%lld,%f,%lld,%Lf\n",
    //     problem->name, PRI(problem->level), PRI(problem->graph->n),PRI(problem->graph->nz),
    //     PRI(problem->graph->max_vtx_cost), gamma, 
    //     PRI(mc_count),TIMER_ELAPSED) ;
    // }
    // fflush(stdout) ;

    
    PERFORMANCE_UPDATE_TIME(MC) ;

    return;
}

VSPFLOAT Reduce_Gamma
(
    VSPFLOAT       gammaIn,
    GR_Reduce_Type redType,
    State           *Z,
    State           *W,
    Options   *options,
    Problem   *problem
)
{
    SET_FUNCTION_COUNTER ;
    VSPFLOAT ret ;

    switch (redType)
    {
        case NO_REDUCE :
            ret = gammaIn ;
            break;

        case USE_MINREDUCE :
            ret = minreduce(Z, W, options, problem) - 1E-3 ;
            break ;

        case REDUCE_SIMPLY :
            ret = gammaIn / 1.3 ;
            break ;

        default:
            #ifdef DEBUG
                ERR("Invalid gamma reduction type") ;
            #endif
            ret = -1 ;
            break ;
    }

    return ret ;
}
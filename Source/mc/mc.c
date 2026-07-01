#include "core.h"
#include "mc_internal.h"
#include "profiler.h"
#include "utils.h"
#include <stdio.h>

#define EPS (1E-10)

void mountain_climbing
(
    State    *X,
    State    *Y,
    Problem  *problem,
    Options  *options,
    VSPFLOAT  gamma,
    bool      is_gamma_org
)
{
    TIMER_START(T_MC) ;

SET_FUNCTION_COUNTER ;
#ifdef DEBUG
        if(is_gamma_org) VSP_ASSERT(equal(gamma,problem->graph->max_vtx_cost)) ;
        else             VSP_ASSERT(gamma < problem->graph->max_vtx_cost) ;

        /* Check input */
        check_state(X, problem);
        check_state(Y, problem);

        VSPFLOAT obj_dbg ;
        obj_dbg = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj_dbg ;
#endif

    /* Local vars */
    VSPFLOAT obj, old_obj, new_obj, f_xstar_y, f_x_ystar ;
    VSPINT mc_count ;
    VSPINT *c = problem->graph->c ;
    VSPINT  n = problem->graph->n ;
    VSPINT  u = problem->u ;
    State *Xstar, *Ystar ;
    char optimal_var;

    /* ====================================================================== */
    // Start mc: either start from lowest cost variable, or start with both
    /* ====================================================================== */
    if(options->start_mc_both_vars)
    {
        obj = compute_obj(X, Y, problem, gamma);
        State *Xstar = Clone_State(X) ;
        State *Ystar = Clone_State(Y) ;

        // Optimize over x to get x_star and f(x*, y)
        f_xstar_y = LP_solver_4
                    (Xstar, Y->Az, c, u, gamma, problem, options, 'x');
        f_xstar_y += Y->weight;

        // Optimize over y to get y_star and f(x, y*)
        f_x_ystar = LP_solver_4
                    (Ystar, X->Az, c, u, gamma, problem, options, 'y');
        f_x_ystar += X->weight;

        /* Adjust objective value to prevent cycling */
        f_xstar_y = MAX(f_xstar_y, obj);
        f_x_ystar = MAX(f_x_ystar, obj);

        if (f_xstar_y > f_x_ystar + f_x_ystar * EPS)
        {
            copy_state(Xstar, X);
            optimal_var = 'x';
        }
        else
        {
            copy_state(Ystar, Y);
            optimal_var = 'y';
        }

        new_obj = (optimal_var == 'x') ? f_xstar_y : f_x_ystar;

        /* Free memory */
        Destroy_State(Xstar);
        Destroy_State(Ystar);
    }
    else
    {
        if(X->weight <= Y->weight)
        {
            // Optimize over x
            new_obj     = LP_solver_4(X, Y->Az, c, u,gamma,problem,options,'x');
            new_obj    += Y->weight;
            optimal_var = 'x';
        }
        else
        {
            // Optimize over y
            new_obj      = LP_solver_4(Y, X->Az, c,u,gamma,problem,options,'y');
            new_obj     += X->weight;
            optimal_var  = 'y';
        }
    }

    /*** Orthogonal Phase *****************************************************/
    old_obj = new_obj ;
    mc_count = 0 ;
    while (true)
    {
        if (optimal_var == 'x')
        {
            /* Optimize over y */
            new_obj = LP_solver_4(Y, X->Az, c, u, gamma, problem, options, 'y');
            new_obj += X->weight ;
            optimal_var = 'y';
        }
        else
        {
            /* Optimize over x */
            new_obj = LP_solver_4(X, Y->Az, c, u, gamma, problem, options, 'x');
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

    TIMER_STOP(T_MC) ;
    return;
}
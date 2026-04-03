#include "shared_api.h"
#include "mc_internal.h"
#include "shared_functions.h"
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
    static long long unsigned count = 0 ;
    count++ ;

    // printf("\nStarting mc with gamma_org=%lld", (long long) is_gamma_org ) ;
    // printf("Value of gamma   : %f", gamma) ;
    // printf("Val of gamma_org : %lld", (long long) problem->Graph->max_vtx_cost) ;
    fflush(stdout) ;

    #ifdef DEBUG
        if(is_gamma_org)
        {
            VSP_ASSERT(equal(gamma, problem->Graph->max_vtx_cost)) ;
        }
        else
        {
            VSP_ASSERT(gamma < problem->Graph->max_vtx_cost) ;
        }
    #endif


    TIMER_START ;

    #ifdef DEBUG
        /* Check input */
        check_state(Y, problem);
        check_state(X, problem);

        VSPFLOAT obj_dbg ;
        obj_dbg = compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    /* Local vars */
    VSPINT *c  = problem->Graph->c;
    VSPINT bal = problem->u;
    VSPINT n   = problem->Graph->n;

    /* Allocate memory */
    State *Xstar = Construct_State(n);
    State *Ystar = Construct_State(n);

    /* Initialize Xstar and Ystar as X and Y */
    copy_state(X, Xstar);
    copy_state(Y, Ystar);

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

            // TEMP
            // VSPFLOAT obj_dbg = compute_obj(State *X, State *Y, Problem *problem, VSPFLOAT gamma)
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
                    c, problem->u, gamma, problem->Graph->total_vtx_cost,
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
                    c, problem->u, gamma, problem->Graph->total_vtx_cost,
                    options
                ) ;
                if(is_opt) break ;
            }
        }

        old_obj = new_obj;
    }
    TIMER_END ;

    // printf("level: %lld  n: %lld  nz: %lld gamma_org: %lld gamma_used: %f mc_count: %lld  mc_time: %Lf\n",
    // PRI(problem->level), PRI(n), PRI(problem->Graph->nz), PRI(problem->Graph->max_vtx_cost),
    // gamma, PRI(mc_count),
    // TIMER_ELAPSED ) ;

    // if(!is_gamma_org)
    // {
    //     printf("%s,%lld,%lld,%lld,%lld,%f,%lld,%Lf\n",
    //     problem->name, PRI(problem->level), PRI(problem->Graph->n),PRI(problem->Graph->nz),
    //     PRI(problem->Graph->max_vtx_cost), gamma, 
    //     PRI(mc_count),TIMER_ELAPSED) ;
    // }
    // fflush(stdout) ;

    
    PERFORMANCE_UPDATE_TIME(MC) ;

    return;
}


/*  Starting from (x, y), we have two choices on where to go:
        Start by optimizing over x
        Start by optimizing over y
    Traditional mc solves one LP in each direction, and picks the best direction
    (or the diagonal point). This version of mc solves "numLPs" LPs before
    making a final decision on the direction. For example, if numLPs = 3, this
    function does
        Start with x: (x0, y0) --> (x1, y0)  --> (x1, y1)   --> (x2, y1)
        Start with y: (x0, y0) --> (x0, y1*) --> (x1*, y1*) --> (x*, y2*)
    and compares (x2, y1) with (x*, y2*) */
// void mountain_climbing_v2
// (
//     State                  *X,
//     State                  *Y,
//     Problem          *problem,
//     Options          *options,
//     VSPFLOAT                gamma,
//     VSPINT                 numLPs /* >= 0. If 0, we launch directly into ortho*/
// )
// {
// #ifdef DEBUG
//     static long long unsigned count = 0 ;
//     count++;

//     /* Check input */
//     check_state(Y, problem);
//     check_state(X, problem);

//     VSPFLOAT obj_dbg ;
//     obj_dbg = compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
//     (void) obj_dbg ;
// #endif
//     VSP_ASSERT(numLPs >= 0) ;

//     /* Local vars */
//     VSPINT *c = problem->Graph->c;
//     VSPINT bal = problem->u;
//     VSPINT n = problem->Graph->n;


//     /* Allocate memory */
//     State *Xstar = Construct_State(n);
//     State *Ystar = Construct_State(n);

//     /* Initialize Xstar and Ystar as X and Y */
//     copy_state(X, Xstar);
//     copy_state(Y, Ystar);

//     VSPINT loop = 0;
//     char optimal_var;

//     /*** Diagonal Phase *******************************************************/
//     /* Save a copy of the starting point */
//     copy_state(X, Xstar) ;
//     copy_state(Y, Ystar) ;

//     /* Start by solving over x */
//     VSPINT LPcounter = 0 ;
//     while( 1 )
//     {

//         LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//         LPcounter++ ;

//         if(LPcounter == numLPs) break ;

//         LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//         LPcounter++ ;

//         if(LPcounter == numLPs) break ;
//     }

//     /* Solve over y */
//     LPcounter = 0 ;
//     while( 1 )
//     {
//         LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//         LPcounter++ ;
//         #ifdef DEBUG
//             UNUSED VSPFLOAT objDbg = compute_obj(X, Y, problem, gamma) ;

//         #endif
//         if(LPcounter == numLPs) break ;

//         LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//         LPcounter++ ;
//         #ifdef DEBUG
//             objDbg = compute_obj(X, Y, problem, gamma) ;
//         #endif
//         if(LPcounter == numLPs) break ;
//     }

//     /* Determine which point is better */
//     VSPFLOAT obj1 = compute_obj(X, Y, problem, gamma) ;
//     VSPFLOAT obj2 = compute_obj(Xstar, Ystar, problem, gamma) ;

//     if(obj2 > obj1)
//     {
//         copy_state(Xstar, X) ;
//         copy_state(Ystar, Y) ;
//     }

//     Destroy_State(Xstar) ;
//     Destroy_State(Ystar) ;

//     /*** Orthogonal Phase *****************************************************/
//     optimal_var = 'x' ;
//     while (true)
//     {
//         loop++;

//         if (optimal_var == 'x')
//         {
//             /* Optimize over y */
//             LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//             optimal_var = 'y';
//         }
//         else
//         {
//             /* Optimize over x */
//             LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//             optimal_var = 'x';
//         }

//         if(optimal_var == 'y')
//         {
//             if( LP_is_optimal(X, Y->Az, gamma, problem, options) )
//             {
//                 break ;
//             }
//         }
//         else
//         {
//             if( LP_is_optimal(Y, X->Az, gamma, problem, options) )
//             {
//                 break ;
//             }
//         }
//     }
//     return;
// }



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
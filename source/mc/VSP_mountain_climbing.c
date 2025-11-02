#include "shared_api.h"
#include "mc_internal.h"
#include "Stopwatches.h"
#include "shared_functions.h"

#define VSP_EPS (1E-10)

void VSP_mountain_climbing
(
    VSP_state               *X,
    VSP_state               *Y,
    VSP_problem       *problem,
    VSP_options       *options,
    VSPFLOAT             gamma,
    bool          is_gamma_org
)
{
    static long long unsigned count = 0 ;
    count++ ;

    printf("Starting mc with gamma_org=%lld\n", (long long) is_gamma_org ) ;
    printf("Value of gamma   : %f", gamma) ;
    printf("Val of gamma_org : %lld", (long long) problem->Graph->max_vtx_cost) ;
    fflush(stdout) ;

    #ifdef VSP_DEBUG
        if(is_gamma_org)
        {
            VSP_ASSERT(VSP_equal(gamma, problem->Graph->max_vtx_cost)) ;
        }
        else
        {
            VSP_ASSERT(gamma < problem->Graph->max_vtx_cost) ;
        }
    #endif

    if( is_gamma_org )
    {
        Stopwatches_Start(problem->Stopwatches, "MC_OrgGamma") ;
    }
    else
    {
        Stopwatches_Start(problem->Stopwatches, "MC_RedGamma") ;
    }

    TIMER_START ;

    #ifdef VSP_DEBUG
        /* Check input */
        VSP_check_state(Y, problem);
        VSP_check_state(X, problem);

        VSPFLOAT obj_dbg ;
        obj_dbg = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    /* Local vars */
    VSPINT *c  = problem->Graph->c;
    VSPINT bal = problem->u;
    VSPINT n   = problem->Graph->n;

    /* Allocate memory */
    VSP_state *Xstar = VSP_construct_state(n);
    VSP_state *Ystar = VSP_construct_state(n);

    /* Initialize Xstar and Ystar as X and Y */
    VSP_copy_state(X, Xstar);
    VSP_copy_state(Y, Ystar);

    VSPFLOAT obj = VSP_compute_obj(X, Y, problem, gamma);
    VSPINT loop = 0;
    char optimal_var;

    /*** Diagonal Phase *******************************************************/
    while (true)
    {
        loop++;

        /* Optimize over x to get x_star */
        VSPFLOAT f_xstar_y = VSP_LP_solver_4(Xstar, Y->Az, c, bal,
                                                 gamma, problem, options, 'x') ;
        f_xstar_y += Y->weight;

        /* Optimize over y to get y_star */
        VSPFLOAT f_x_ystar = VSP_LP_solver_4(Ystar, X->Az, c, bal,
                                                 gamma, problem, options, 'y') ;
        f_x_ystar += X->weight;

        /* Adjust objective value to prevent cycling */
        f_xstar_y = MAX(f_xstar_y, obj);
        f_x_ystar = MAX(f_x_ystar, obj);

        /* Compare f(x*, y*) vs f(x_star, y) vs f(x, y_star) */
        VSPFLOAT f_xstar_ystar;
        f_xstar_ystar = VSP_compute_obj(Xstar, Ystar, problem, gamma);

        if (f_xstar_ystar > MAX(f_xstar_y, f_x_ystar) + MAX(f_xstar_y, f_x_ystar) * VSP_EPS)
        {
            /* Set (x, y) <- (x_star, y_star) and repeat */
            /* TODO speed: Don't have 2 states and copy. Just maintain 4 states
               throughout an entire level                                     */
            VSP_copy_state(Xstar, X);
            VSP_copy_state(Ystar, Y);

            obj = f_xstar_ystar;
            continue;
        }
        else if (f_xstar_y > f_x_ystar + f_x_ystar * VSP_EPS)
        {
            VSP_copy_state(Xstar, X);
            optimal_var = 'x';
            break ;
        }
        else
        {
            VSP_copy_state(Ystar, Y);
            optimal_var = 'y';
            break;
        }
    }

    /* Free memory */
    VSP_destruct_state(Xstar);
    VSP_destruct_state(Ystar);

    /*** Orthogonal Phase *****************************************************/
    VSPFLOAT old_obj = VSP_compute_obj(X, Y, problem, gamma);
    while (true)
    {
        VSPFLOAT new_obj;
        loop++;

        if (optimal_var == 'x')
        {
            /* Optimize over y */
            new_obj = VSP_LP_solver_4(Y, X->Az, c, bal, gamma, problem, options, 'y') ;
            new_obj += X->weight;
            optimal_var = 'y';
        }
        else
        {
            /* Optimize over x */
            new_obj = VSP_LP_solver_4(X, Y->Az, c, bal, gamma, problem, options, 'x') ;
            new_obj += Y->weight;
            optimal_var = 'x';

            // TEMP
            // VSPFLOAT obj_dbg = VSP_compute_obj(VSP_state *X, VSP_state *Y, VSP_problem *problem, VSPFLOAT gamma)
        }

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
    PERFORMANCE_UPDATE_TIME(MC) ;

    if( is_gamma_org )
    {
        Stopwatches_Stop(problem->Stopwatches, "MC_OrgGamma") ;
    }
    else
    {
        Stopwatches_Stop(problem->Stopwatches, "MC_RedGamma") ;
    }
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
// void VSP_mountain_climbing_v2
// (
//     VSP_state                  *X,
//     VSP_state                  *Y,
//     VSP_problem          *problem,
//     VSP_options          *options,
//     VSPFLOAT                gamma,
//     VSPINT                 numLPs /* >= 0. If 0, we launch directly into ortho*/
// )
// {
// #ifdef VSP_DEBUG
//     static long long unsigned count = 0 ;
//     count++;

//     /* Check input */
//     VSP_check_state(Y, problem);
//     VSP_check_state(X, problem);

//     VSPFLOAT obj_dbg ;
//     obj_dbg = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
//     (void) obj_dbg ;
// #endif
//     VSP_ASSERT(numLPs >= 0) ;

//     /* Local vars */
//     VSPINT *c = problem->Graph->c;
//     VSPINT bal = problem->u;
//     VSPINT n = problem->Graph->n;


//     /* Allocate memory */
//     VSP_state *Xstar = VSP_construct_state(n);
//     VSP_state *Ystar = VSP_construct_state(n);

//     /* Initialize Xstar and Ystar as X and Y */
//     VSP_copy_state(X, Xstar);
//     VSP_copy_state(Y, Ystar);

//     VSPINT loop = 0;
//     char optimal_var;

//     /*** Diagonal Phase *******************************************************/
//     /* Save a copy of the starting point */
//     VSP_copy_state(X, Xstar) ;
//     VSP_copy_state(Y, Ystar) ;

//     /* Start by solving over x */
//     VSPINT LPcounter = 0 ;
//     while( 1 )
//     {

//         VSP_LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//         LPcounter++ ;

//         if(LPcounter == numLPs) break ;

//         VSP_LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//         LPcounter++ ;

//         if(LPcounter == numLPs) break ;
//     }

//     /* Solve over y */
//     LPcounter = 0 ;
//     while( 1 )
//     {
//         VSP_LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//         LPcounter++ ;
//         #ifdef VSP_DEBUG
//             UNUSED VSPFLOAT objDbg = VSP_compute_obj(X, Y, problem, gamma) ;

//         #endif
//         if(LPcounter == numLPs) break ;

//         VSP_LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//         LPcounter++ ;
//         #ifdef VSP_DEBUG
//             objDbg = VSP_compute_obj(X, Y, problem, gamma) ;
//         #endif
//         if(LPcounter == numLPs) break ;
//     }

//     /* Determine which point is better */
//     VSPFLOAT obj1 = VSP_compute_obj(X, Y, problem, gamma) ;
//     VSPFLOAT obj2 = VSP_compute_obj(Xstar, Ystar, problem, gamma) ;

//     if(obj2 > obj1)
//     {
//         VSP_copy_state(Xstar, X) ;
//         VSP_copy_state(Ystar, Y) ;
//     }

//     VSP_destruct_state(Xstar) ;
//     VSP_destruct_state(Ystar) ;

//     /*** Orthogonal Phase *****************************************************/
//     optimal_var = 'x' ;
//     while (true)
//     {
//         loop++;

//         if (optimal_var == 'x')
//         {
//             /* Optimize over y */
//             VSP_LP_solver_4(Y, X->Az, c, 1, bal, gamma, problem, options, 'y');
//             optimal_var = 'y';
//         }
//         else
//         {
//             /* Optimize over x */
//             VSP_LP_solver_4(X, Y->Az, c, 1, bal, gamma, problem, options, 'x');
//             optimal_var = 'x';
//         }

//         if(optimal_var == 'y')
//         {
//             if( VSP_LP_is_optimal(X, Y->Az, gamma, problem, options) )
//             {
//                 break ;
//             }
//         }
//         else
//         {
//             if( VSP_LP_is_optimal(Y, X->Az, gamma, problem, options) )
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
    VSP_state           *Z,
    VSP_state           *W,
    VSP_options   *options,
    VSP_problem   *problem
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
            ret = VSP_minreduce(Z, W, options, problem) - 1E-3 ;
            break ;

        case REDUCE_SIMPLY :
            ret = gammaIn / 1.3 ;
            break ;

        default:
            #ifdef VSP_DEBUG
                ERR("Invalid gamma reduction type") ;
            #endif
            ret = -1 ;
            break ;
    }

    return ret ;
}
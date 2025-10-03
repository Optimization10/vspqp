#include "Stopwatches.h"
#include "shared_api.h"
#include "ds_api.h"
#include "gr_api.h"
#include "flips_api.h"
#include "FM_api.h"
#include "FM_flips_api.h"
#include "mc_api.h"
#include "flips_public_v4.h"
#include "flips_public_v5.h"
#include "flips_public_v7.h"
#include "string.h"
#include "gr_flips.h"
#include "logger.h"


static void Best_Of_MC_And_Greedy_Flips
(
    VSP_state                *X,
    VSP_state                *Y,
    VSP_problem        *problem,
    VSP_options        *options,
    VSPINT         numFlipsRuns
) ;


void FM_Flips
(
    VSP_state         *X,
    VSP_state         *Y,
    VSP_problem *problem
) ; 

static void Determine_Num_GR
(
    VSPINT             n,
    VSP_options *options,
    VSPINT        *numGR  // output
)
{
    *numGR = 0 ;
    if(options->use_gamma_reduction)
    {
        if(n <= options->max_num_vert_small_graph)
        {
            *numGR = options->num_gr_for_small_graphs ;
        }
        else
        {
            *numGR = options->num_gr_for_large_graphs ;
        }
    }
    return ;
}

// whereToUseFMFlips:
// 0 - FMv2 disabled (or anything besides 1, 2, 3)
// 1 - FMv2 only during first optimization
// 2 - FMv2 during gamma reduction
// 3 - FMv2 in both locations
void Optimize
(
    VSP_state                     *X,
    VSP_state                     *Y,
    VSP_problem             *problem,
    VSP_options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
)
{
    static long long unsigned count = 0 ;
    count++ ;
    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif

    /* 0.1 Local vars */
    VSPINT i, num_gamma_reductions, n = problem->Graph->n, gamma_org = problem->Graph->max_vtx_cost ;

    VSPFLOAT f_x_y, f_xstar_ystar ;

    VSP_state *Xstar = VSP_construct_state( n ) ;
    VSP_state *Ystar = VSP_construct_state( n ) ;

    /* 2. Solve the problem using the original gamma */
    VSPFLOAT gamma = problem->Graph->max_vtx_cost ;


    /*--- Solve VSP ----------------------------------------------------------*/
    /* First refinement */
    LoggerStart(X, Y, problem, "MC And Greedy") ;

    /*------------------------------------------------------------------------*/
    // 1. mountain climbing vs Flips
    Best_Of_MC_And_Greedy_Flips(X, Y, problem, options, 5) ;

    /*------------------------------------------------------------------------*/

    LoggerStop(X, Y, problem) ;

    LoggerStart(X, Y, problem, "Flips") ;
    Flips(X, Y, problem, 9, VERSION_DEFAULT) ;

    LoggerStop(X, Y, problem) ;

    VSP_relaxed_to_vertex_separator(X, Y, problem) ;

    LoggerStart(X, Y, problem, "Metis") ;
    VSPINT n_passes ;
    if(problem->level >= problem->coarsest_level_num / 2) n_passes = 3 ;
    else n_passes = 2 ;

    
    FM_Metis(X, Y, problem, options, n_passes) ;

    LoggerStop(X, Y, problem) ;

    /*------------------------------------------------------------------------*/
    if(whereToUseFMFlips == 1 || whereToUseFMFlips == 3)
    {       
        LoggerStart(X, Y, problem, "FM Flips") ;
        FM_Flips(X, Y, problem) ;
        LoggerStop(X, Y, problem) ;
    }
    /*------------------------------------------------------------------------*/

    Force_To_Binary_Point(X, Y, problem) ;


    #ifdef VSP_DEBUG
        VSPFLOAT obj_dbg ;
        obj_dbg = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    /* Gamma reductions */
    VSPINT loop = 0 ;
    Determine_Num_GR(n, options, &num_gamma_reductions) ;
    VSPINT k = 0 ;

    for(i = 0 ; i < num_gamma_reductions ; i++)
    {
        /*--------------------------------------------------------------------*/
        /* Reduce gamma                                                       */
        /*--------------------------------------------------------------------*/
        if(k == 0)
        {
            gamma = VSP_minreduce(X, Y, options, problem) - 1E-3 ;
            if(gamma < 0) break ;
        }
        else gamma = gamma / 1.3 ;

        /*--------------------------------------------------------------------*/
        /* (x*, y*) <- Solution starting from (x, y), using the reduced gamma */
        /*--------------------------------------------------------------------*/
        VSP_copy_state(X, Xstar) ;
        VSP_copy_state(Y, Ystar) ;

        VSP_mountain_climbing(Xstar, Ystar, problem, options, gamma) ;

        /*--------------------------------------------------------------------*/
        /* (x*, y*) <- Solution starting from (x*, y*), using original gamma  */
        /*--------------------------------------------------------------------*/
        VSP_mountain_climbing(Xstar, Ystar, problem, options, gamma_org) ;
        Flips(Xstar, Ystar, problem, 9, VERSION_DEFAULT) ;

        VSP_relaxed_to_vertex_separator(Xstar, Ystar, problem) ;

        VSPINT n_passes ;
        if(problem->level >= problem->coarsest_level_num / 2) n_passes = 3 ;
        else n_passes = 2 ;

        Stopwatches_Start(problem->Stopwatches, "FM_METIS") ;
        FM_Metis(Xstar, Ystar, problem, options, n_passes) ;
        Stopwatches_Stop(problem->Stopwatches, "FM_METIS") ;

        /*--------------------------------------------------------------------*/
        if(whereToUseFMFlips == 2 || whereToUseFMFlips == 3)
        {          
            FM_Flips(Xstar, Ystar, problem) ;
        }
        /*--------------------------------------------------------------------*/

        /*--------------------------------------------------------------------*/
        /* Swap (x*, y*) with (x, y) if it is beneficial to do so             */
        /*--------------------------------------------------------------------*/
        f_x_y         = VSP_compute_obj(X, Y, problem, gamma_org) ;
        f_xstar_ystar = VSP_compute_obj(Xstar, Ystar, problem, gamma_org  );

        if(f_xstar_ystar > f_x_y + fabs(f_x_y) * 1E-13 )
        {
            Swap_Mem(X, Xstar, sizeof(VSP_state) ) ;
            Swap_Mem(Y, Ystar, sizeof(VSP_state) ) ;
            k = 0 ; // Set k so that we call minreduce again
        }
        else k++ ;
    }
    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif
    /*------------------------------------------------------------------------*/
    // Do a last round of flips to ensure that we have a bona fide separator.
    // The way a separator might not happen: in minreduce, we do a MC which
    // can produce penalties. This is the last thing done do the x, y point
    /*------------------------------------------------------------------------*/
    LoggerStart(X, Y, problem, "Final Flips") ;
    Flips(X, Y, problem, 5, VERSION_DEFAULT) ;
    LoggerStop(X, Y, problem) ;

    /*------------------------------------------------------------------------*/
    /* Cleanup                                                                */
    /*------------------------------------------------------------------------*/
    VSP_destruct_state(Xstar) ;
    VSP_destruct_state(Ystar) ;

    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif

    return ;
}

extern struct timespec G_START_TIME, G_END_TIME ;
extern long double G_FLIPS_TIME ;

// Takes in a feasible point for the continuous problem, forces it to a binary
// point, and optimizes it
// Returns a feasible point for the binary problem with objective at least
// as large as it started with after projecting the fractional entry to 0 or 1
void Flips
(
    VSP_state               *X,
    VSP_state               *Y,
    VSP_problem       *problem,
    VSPINT        maxNumPasses,
    VersionNum         version
)
{
    clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;

    SET_FUNCTION_COUNTER ;
    Stopwatches_Start(problem->Stopwatches, "Flips") ;
    Force_To_Binary_Point(X, Y, problem) ;

    VSPINT *x  = X->point->int_part ;
    VSPINT *Ax = X->Az->int_part ;
    VSPINT *y  = Y->point->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT *wgtx = &(X->weight) ;
    VSPINT *wgty = &(Y->weight) ;


    switch (version)
    {
        case VERSION_DEFAULT:
            if(problem->level <= 2)
            {
                Flips_v5(x, Ax, wgtx, y, Ay, wgty, problem, maxNumPasses) ;
            }
            else
            {
                Flips_v4(x, Ax, wgtx, y, Ay, wgty, problem, maxNumPasses) ;
            }
            break ;
        case VERSION_4:
            Flips_v4(x, Ax, wgtx, y, Ay, wgty, problem, maxNumPasses) ;
            break ;
        case VERSION_5:
            Flips_v5(x, Ax, wgtx, y, Ay, wgty, problem, maxNumPasses) ;
            break ;
        default:
            printf("\nError: Did not select correct flips version") ;
            break ;
    }

    Stopwatches_Stop(problem->Stopwatches, "Flips") ;

    clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
    G_FLIPS_TIME += VSP_compute_elapsed_time(G_START_TIME, G_END_TIME) ;
}

void FM_Flips
(
    VSP_state         *X,
    VSP_state         *Y,
    VSP_problem *problem
)
{
    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif

    if(problem->level >= 0)
    {
        // FM_Flips_API_v2(X, Y, problem, options) ; 
        FM_Flips_v7(X, Y, problem, 8) ;
        // FM_Flips_v8(X, Y, problem, 4) ;
        // FM_Flips_v9(X, Y, problem, 4) ;
    }

    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif
}

static void Best_Of_MC_And_Greedy_Flips
(
    VSP_state                *X,
    VSP_state                *Y,
    VSP_problem        *problem,
    VSP_options        *options,
    VSPINT         numFlipsRuns
)
{
    static long long unsigned count = 0 ;
    count++ ;
    
    VSPINT n     = problem->Graph->n ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;


    /*------------------------------------------------------------------------*/
    // Make a copy of the input data
    /*------------------------------------------------------------------------*/
    VSPINT *xcpy, *Axcpy, *ycpy, *Aycpy, wgtxcpy, wgtycpy ;
    xcpy  = malloc(n * sizeof(VSPINT)) ;
    Axcpy = malloc(n * sizeof(VSPINT)) ;
    ycpy  = malloc(n * sizeof(VSPINT)) ;
    Aycpy = malloc(n * sizeof(VSPINT)) ;
    Copy_Data
    (
        X->point->int_part, X->Az->int_part, &(X->weight),
        Y->point->int_part, Y->Az->int_part, &(Y->weight),
        xcpy, Axcpy, &wgtxcpy,
        ycpy, Aycpy, &wgtycpy,
        problem->Graph->n
    ) ;

    // Compute c*x and c*y locally
    wgtxcpy = 0 ;
    wgtycpy = 0 ;
    VSPINT *c = problem->Graph->c ;
    for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
    {
        wgtxcpy += xcpy[i] * c[i] ;
        wgtycpy += ycpy[i] * c[i] ;
    }

    bool abortManyFlipsPasses = false ;
    if(wgtxcpy == 0 || wgtycpy == 0)
    {
        abortManyFlipsPasses = true ;       
    }

    /*------------------------------------------------------------------------*/
    // Do greedy flips, all from the same starting point 
    /*------------------------------------------------------------------------*/

    if(!abortManyFlipsPasses)
    {
        TryFlipsManyStartingPoints
        (
            &xcpy, &Axcpy, &wgtxcpy,
            &ycpy, &Aycpy, &wgtycpy,
            problem, options,
            numFlipsRuns, 9
        ) ;
    }

    /*------------------------------------------------------------------------*/
    // Do Mountain CLimbing in-place in (x, y)
    /*------------------------------------------------------------------------*/
    
    VSP_mountain_climbing(X, Y, problem, options,  gamma) ;

    Force_To_Binary_Point(X, Y, problem) ;

    /*------------------------------------------------------------------------*/
    // Compare MC point with Flips point and choose the best
    /*------------------------------------------------------------------------*/
    VSPINT objAfterMc ;
    Compute_ObjValSimple
    (
        X->point->int_part, X->Az->int_part, X->weight,
        Y->point->int_part, Y->Az->int_part, Y->weight,
        problem, &objAfterMc
    ) ;

    VSPINT objAfterManyFlips ;
    Compute_ObjValSimple
    (
        xcpy, Axcpy, wgtxcpy,
        ycpy, Aycpy, wgtycpy,
        problem, &objAfterManyFlips
    ) ;

    bool objFlipsBetter = (objAfterManyFlips > objAfterMc) ;
    bool objFlipsSame   = (objAfterManyFlips == objAfterMc) ;
    bool balFlipsBeter  = (IMBALANCE(wgtxcpy, wgtycpy) < IMBALANCE( X->weight, Y->weight) );
    if( !abortManyFlipsPasses &&(objFlipsBetter || (objFlipsSame && balFlipsBeter) ) )
    {
        // Overwrite (x, y) with the point from flips
        Copy_Data
        (
            xcpy, Axcpy, &wgtxcpy,
            ycpy, Aycpy, &wgtycpy,
            X->point->int_part, X->Az->int_part, &(X->weight),
            Y->point->int_part, Y->Az->int_part, &(Y->weight),
            n
        ) ;
        X->point->ind_frac = EMPTY ;
        X->Az->ind_frac_z  = EMPTY ;
        Y->point->ind_frac = EMPTY ;
        Y->Az->ind_frac_z  = EMPTY ;

        #ifdef VSP_DEBUG
            VSP_check_states(X, Y, problem) ;
            VSPINT obj = (VSPINT) VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
            if(obj != objAfterManyFlips )
            {
                ERR("PROBLEM") ;
            }
        #endif
    }

    /*------------------------------------------------------------------------*/
    // Cleanup and return
    /*------------------------------------------------------------------------*/

    free(xcpy)  ;
    free(Axcpy) ;
    free(ycpy)  ;
    free(Aycpy) ;
}
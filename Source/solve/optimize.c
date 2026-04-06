#include "shared_api.h"
#include "ds_api.h"
#include "gr_api.h"
#include "flips_api.h"
#include "mc_api.h"
#include "FM_api.h"
#include "flips_public_v4.h"
#include "shared_functions.h"
#include "logger.h"
#include <stdio.h>
#include "solve_internal.h"


void TryFlipsManyStartingPoints
(
    VSPINT                 **xBestHome,
    VSPINT                **AxBestHome,
    VSPINT               *wgtxBestHome,
    VSPINT                 **yBestHome,
    VSPINT                **AyBestHome,
    VSPINT               *wgtyBestHome,
    Problem               *problem,
    Options               *options,
    VSPINT                 numFlipRuns,
    VSPINT            numPassesInFlips
) ;

static void Best_Of_MC_And_Greedy_Flips
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options,
    VSPINT  numFlipsRuns
)
{
    static long long unsigned count = 0 ;
    count++ ;
    
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;


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
        problem->graph->n
    ) ;

    // Compute c*x and c*y locally
    VSPINT *c = problem->graph->c ;


    bool abortManyFlipsPasses = false ;
    wgtxcpy = State_Get_Wgt_Int_Part(X, c) ;
    wgtycpy = State_Get_Wgt_Int_Part(Y, c) ;
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
    
    mountain_climbing(X, Y, problem, options,  gamma, true) ;

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

        #ifdef DEBUG
            check_states(X, Y, problem) ;
            VSPINT obj = (VSPINT) compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
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

static void Determine_Num_GR
(
    VSPINT             n,
    Options *options,
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


// Single an early exit if f(x, y) >= c(V) - alpha
// The logic is that if (X, Y) represent a vertex separator
// partition, and the separator only contains a node of minimal weight, then
// X, Y is globally optimal, or nearly globally optimal
// Note that we are not taking into account the penalty term
bool Early_Exit(State *X, State *Y, Problem *P)
{
    VSPFLOAT obj = compute_obj(X, Y, P, P->graph->max_vtx_cost) ;
    bool ret = false ;
    if(lteq(P->graph->total_vtx_cost - P->graph->min_vtx_cost, obj))
    {
        ret = true ;
    }
    return ret ;
}


void Optimize
(
    State                     *X,
    State                     *Y,
    Problem             *problem,
    Options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
)
{
    static long long unsigned count = 0 ;
    count++ ;
    #ifdef DEBUG
        check_states(X, Y, problem) ;
    #endif

    /* 0.1 Local vars */
    VSPINT i, k, num_gamma_reductions, n = problem->graph->n,
    gamma_org = problem->graph->max_vtx_cost ;

    VSPFLOAT f_x_y, f_xstar_ystar ;

    State *Xstar = Construct_State( n ) ;
    State *Ystar = Construct_State( n ) ;

    /* 2. Solve the problem using the original gamma */
    VSPFLOAT gamma = problem->graph->max_vtx_cost ;

    /*--- Solve VSP ----------------------------------------------------------*/
    // 1. mountain climbing vs Flips
    Best_Of_MC_And_Greedy_Flips(X, Y, problem, options, 5) ;

    if(Early_Exit(X, Y, problem))
    {
        goto EARLY_EXIT ;
    }

    // 2. Flips
    LoggerStart(X, Y, problem, "Flips") ;
    Flips(X, Y, problem, 9) ;
    LoggerStop(X, Y, problem) ;

    // 3. Get vertex separator in preparation for FM_METIS
    relaxed_to_vertex_separator(X, Y, problem) ;
    if(Early_Exit(X, Y, problem))
    {
        goto EARLY_EXIT ;
    }

    // FM METIS
    LoggerStart(X, Y, problem, "Metis") ;
    VSPINT n_passes ;
    if(problem->level >= problem->coarsest_level_num / 2) n_passes = 3 ;
    else n_passes = 2 ;
    FM_Metis(X, Y, problem, options, n_passes) ;
    LoggerStop(X, Y, problem) ;

    if(Early_Exit(X, Y, problem))
    {
        goto EARLY_EXIT ;
    }
    /*------------------------------------------------------------------------*/


    #ifdef DEBUG
        VSPFLOAT obj_dbg ;
        obj_dbg = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    /* Gamma reductions */
    k = 0 ;
    Determine_Num_GR(n, options, &num_gamma_reductions) ;
    for(i = 0 ; i < num_gamma_reductions ; i++)
    {
        /*--------------------------------------------------------------------*/
        /* Reduce gamma                                                       */
        /*--------------------------------------------------------------------*/
        if(k == 0)
        {
            gamma = minreduce(X, Y, options, problem) - 1E-3 ;
            if(gamma < 0) break ;
        }
        else gamma = gamma / 1.3 ;

        /*--------------------------------------------------------------------*/
        /* (x*, y*) <- Solution starting from (x, y), using the reduced gamma */
        /*--------------------------------------------------------------------*/
        copy_state(X, Xstar) ;
        copy_state(Y, Ystar) ;

        mountain_climbing(Xstar, Ystar, problem, options, gamma, false) ;

        /*--------------------------------------------------------------------*/
        /* (x*, y*) <- Solution starting from (x*, y*), using original gamma  */
        /*--------------------------------------------------------------------*/
        mountain_climbing(Xstar, Ystar, problem, options, gamma_org, true) ;
        Flips(Xstar, Ystar, problem, 9) ;

        relaxed_to_vertex_separator(Xstar, Ystar, problem) ;

        VSPINT n_passes ;
        if(problem->level >= problem->coarsest_level_num / 2) n_passes = 3 ;
        else n_passes = 2 ;

        FM_Metis(Xstar, Ystar, problem, options, n_passes) ;



        /*--------------------------------------------------------------------*/
        /* Swap (x*, y*) with (x, y) if it is beneficial to do so             */
        /*--------------------------------------------------------------------*/
        f_x_y         = compute_obj(X, Y, problem, gamma_org) ;
        f_xstar_ystar = compute_obj(Xstar, Ystar, problem, gamma_org  );

        if(f_xstar_ystar > f_x_y + fabs(f_x_y) * 1E-13 )
        {
            Swap_Mem(X, Xstar, sizeof(State) ) ;
            Swap_Mem(Y, Ystar, sizeof(State) ) ;
            k = 0 ; // Set k so that we call minreduce again

            if(Early_Exit(X, Y, problem))
            {
                goto EARLY_EXIT ;
            }
        }
        else k++ ;
    }
    #ifdef DEBUG
        check_states(X, Y, problem) ;
    #endif
    /*------------------------------------------------------------------------*/
    // Do a last round of flips to ensure that we have a bona fide separator.
    // The way a separator might not happen: in minreduce, we do a MC which
    // can produce penalties. This is the last thing done do the x, y point
    /*------------------------------------------------------------------------*/
    LoggerStart(X, Y, problem, "Final Flips") ;
    Flips(X, Y, problem, 5) ;
    LoggerStop(X, Y, problem) ;

    /*------------------------------------------------------------------------*/
    /* Cleanup                                                                */
    /*------------------------------------------------------------------------*/
    EARLY_EXIT:

    Destroy_State(Xstar) ;
    Destroy_State(Ystar) ;

    #ifdef DEBUG
        check_states(X, Y, problem) ;
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
    State               *X,
    State               *Y,
    Problem       *problem,
    VSPINT        maxNumPasses
)
{
    clock_gettime(CLOCK_REALTIME, &G_START_TIME) ;

    SET_FUNCTION_COUNTER ;
    Force_To_Binary_Point(X, Y, problem) ;

    VSPINT *x  = X->point->int_part ;
    VSPINT *Ax = X->Az->int_part ;
    VSPINT *y  = Y->point->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT *wgtx = &(X->weight) ;
    VSPINT *wgty = &(Y->weight) ;

    bool deep_search = (problem->level <= 2) ;
    Flips_v999(x, Ax, wgtx, y, Ay, wgty, problem, maxNumPasses, deep_search) ;


    clock_gettime(CLOCK_REALTIME, &G_END_TIME) ;
    G_FLIPS_TIME += compute_elapsed_time(G_START_TIME, G_END_TIME) ;
}


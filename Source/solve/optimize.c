#include "solve_internal.h"
#include "mc_api.h"
#include "FM_api.h"
#include "greedy_flips_api.h"
#include "utils.h"




void Best_Of_MC_And_Greedy_Flips
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
        Multi_Start_Flips
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

    Round_Down_Fractions(X, Y, problem) ;

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
static bool Early_Exit(State *X, State *Y, Problem *P)
{
    VSPFLOAT obj = compute_obj(X, Y, P, P->graph->max_vtx_cost) ;
    bool ret = false ;
    if(lteq(P->graph->total_vtx_cost - P->graph->min_vtx_cost, obj))
    {
        ret = true ;
    }
    return ret ;
}


static void Gamma_Reduction_Loop
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       gamma_org
)
{
    VSPINT i, k = 0, num_gr;
    VSPFLOAT gamma, f_x_y, f_xstar_ystar;

    Determine_Num_GR(problem->graph->n, options, &num_gr);

    for (i = 0; i < num_gr; i++)
    {
        if (k == 0)
        {
            gamma = minreduce(X, Y, options, problem) - 1E-3;
            if (gamma < 0) break;
        }
        else gamma = gamma / 1.3;

        State *Xc = Clone_State(X);
        State *Yc = Clone_State(Y);

        mountain_climbing(Xc, Yc, problem, options, gamma, false);
        mountain_climbing(Xc, Yc, problem, options, gamma_org, true);
        Greedy_Flips(Xc, Yc, problem, 9, (problem->level <= 2));
        relaxed_to_vertex_separator(Xc, Yc, problem);

        VSPINT np;
        if (problem->level >= problem->coarsest_level_num / 2) np = 3;
        else np = 2;
        FM_Metis(Xc, Yc, problem, options, np);

        f_x_y         = compute_obj(X,  Y,  problem, gamma_org);
        f_xstar_ystar = compute_obj(Xc, Yc, problem, gamma_org);

        if (f_xstar_ystar > f_x_y + fabs(f_x_y) * 1E-13)
        {
            Swap_Mem(X, Xc, sizeof(State));
            Swap_Mem(Y, Yc, sizeof(State));
            k = 0;

            if (Early_Exit(X, Y, problem))
            {
                Destroy_State(Xc);
                Destroy_State(Yc);
                return;
            }
        }
        else k++;

        Destroy_State(Xc);
        Destroy_State(Yc);
    }
}


void Optimize
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options
)
{
    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    VSPINT gamma_org = problem->graph->max_vtx_cost;

    // Phase 1: MC vs Flips — run both on same input, keep better
    Best_Of_MC_And_Greedy_Flips(X, Y, problem, options, 5);

    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 2: Greedy flips
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2));

    // Phase 3: FM
    relaxed_to_vertex_separator(X, Y, problem);
    VSPINT n_passes;
    if (problem->level >= problem->coarsest_level_num / 2) n_passes = 3;
    else n_passes = 2;
    FM_Metis(X, Y, problem, options, n_passes);

    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 4: Gamma reductions
    Gamma_Reduction_Loop(X, Y, problem, options, gamma_org);

    // Phase 5: Final cleanup flips
    Greedy_Flips(X, Y, problem, 5, (problem->level <= 2));

    CLEANUP:

    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    return;
}
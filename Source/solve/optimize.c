#include "solve_internal.h"
#include "mc.h"
#include "FM.h"
#include "greedy_flips.h"
#include "utils.h"




void Best_Of_MC_And_Greedy_Flips
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       numFlipsRuns
)
{
    VSPINT gamma = problem->graph->max_vtx_cost;
    VSPINT *c = problem->graph->c;

    // Clone for flips path
    State *Xf = Clone_State(X);
    State *Yf = Clone_State(Y);

    VSPINT wgtx = State_Get_Wgt_Int_Part(Xf, c);
    VSPINT wgty = State_Get_Wgt_Int_Part(Yf, c);

    // Flips path
    if (wgtx != 0 && wgty != 0)
    {
        Multi_Start_Flips(Xf, Yf, problem, options, numFlipsRuns, 9) ;
    }

    // MC path (in-place on X, Y)
    mountain_climbing(X, Y, problem, options, gamma, true);


    // Compare
    VSPFLOAT objMC, objFlips ;
    objMC    = compute_obj(X, Y   , problem, (VSPFLOAT) gamma) ;
    objFlips = compute_obj(Xf, Yf , problem, (VSPFLOAT) gamma) ;

    bool flipsBetter    = (objFlips > objMC);
    bool flipsSame      = (objFlips == objMC);
    bool flipsBetterBal = (IMBALANCE(Xf->weight, Yf->weight) < IMBALANCE(X->weight, Y->weight));

    if (wgtx != 0 && wgty != 0 && (flipsBetter || (flipsSame && flipsBetterBal)))
    {
        Swap_Mem(X, Xf, sizeof(State));
        Swap_Mem(Y, Yf, sizeof(State));
    }

    Destroy_State(Xf);
    Destroy_State(Yf);
}

static void Determine_Num_GR
(
    VSPINT             n,
    Options *options,
    VSPINT        *numGR  // output
)
{
    *numGR = 0 ;
    if(n <= options->max_num_vert_small_graph)
    {
        *numGR = options->num_gr_for_small_graphs ;
    }
    else
    {
        *numGR = options->num_gr_for_large_graphs ;
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
        Greedy_Flips(Xc, Yc, problem, 9, (problem->level <= 2), false);

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
    static int count = 0 ;
    count++ ;

    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    VSPINT gamma_org = problem->graph->max_vtx_cost;

    // Phase 1: MC vs Flips — run both on same input, keep better
    // Phase 1: MC
    mountain_climbing(X, Y, problem, options, problem->graph->max_vtx_cost, true);
    

    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 2: Greedy flips
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2), false);

    // Phase 3: FM
    VSPINT n_passes;
    if (problem->level >= problem->coarsest_level_num / 2) n_passes = 3;
    else n_passes = 2;
    FM_Metis(X, Y, problem, options, n_passes);

    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 4: MC
    mountain_climbing
    (
        X, Y, problem, options, problem->graph->max_vtx_cost, true
    ) ;

    // Phase 5: Gamma reductions
    Gamma_Reduction_Loop(X, Y, problem, options, gamma_org);

    // Phase 6: Final cleanup flips
    Greedy_Flips(X, Y, problem, 5, true, true);

    CLEANUP:

    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    return;
}
#include "solve_internal.h"
#include "utils.h"
#include "greedy_flips_api.h"
#include "FM_api.h"

// Run a strategy on a cloned point. If it beats the original, swap.
// Returns true if the clone was better.
bool Try_And_Keep_If_Better
(
    State   *X,
    State   *Y,
    Problem *problem,
    void   (*strategy)(State*, State*, Problem*, Options*),
    Options *options
)
{
    State *Xc = Clone_State(X);
    State *Yc = Clone_State(Y);
    VSPINT gamma = problem->graph->max_vtx_cost;

    strategy(Xc, Yc, problem, options);

    VSPFLOAT obj_orig = compute_obj(X,  Y,  problem, gamma);
    VSPFLOAT obj_new  = compute_obj(Xc, Yc, problem, gamma);

    bool better = (obj_new > obj_orig);
    if (better)
    {
        Swap_Mem(X, Xc, sizeof(State));
        Swap_Mem(Y, Yc, sizeof(State));
    }

    Destroy_State(Xc);
    Destroy_State(Yc);
    return better;
}

void strategy_mc_and_flips
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    Best_Of_MC_And_Greedy_Flips(X, Y, problem, options, 5);
}

void strategy_flips
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2));
}

void strategy_fm
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    relaxed_to_vertex_separator(X, Y, problem);
    VSPINT n_passes;
    if (problem->level >= problem->coarsest_level_num / 2) n_passes = 3;
    else n_passes = 2;
    FM_Metis(X, Y, problem, options, n_passes);
}

void strategy_centroid_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options, true, 0);
    Round_Down_Fractions(X, Y, problem);
}
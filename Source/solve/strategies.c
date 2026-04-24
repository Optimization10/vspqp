#include "VSPQP.h"
#include "solve_internal.h"
#include "utils.h"
#include "greedy_flips.h"
#include "FM.h"

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

void strategy_centroid_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = centroid_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}

void strategy_bisection_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = bisection_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}

void strategy_random_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = random_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}
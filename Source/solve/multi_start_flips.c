
#include "solve_internal.h"
#include "greedy_flips.h"
#include "utils.h"


void Multi_Start_Flips
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       numFlipRuns,
    VSPINT       numPassesInFlips
)
{
    VSPINT n = problem->graph->n;
    VSPINT gamma = problem->graph->max_vtx_cost;

    VSPFLOAT objBest ;
    VSPFLOAT objStar ;
    objBest = compute_obj(X, Y, problem, (VSPFLOAT) gamma) ;

    VSPINT imbBest = IMBALANCE(X->weight, Y->weight);

    State *Xstar = NULL;
    State *Ystar = NULL;
    State *Xbest = NULL;
    State *Ybest = NULL;

    bool foundBetter = false;

    for (VSPINT p = 0; p < numFlipRuns; p++)
    {
        Xstar = Clone_State(X);
        Ystar = Clone_State(Y);

        Greedy_Flips(Xstar, Ystar, problem, numPassesInFlips, true, true);

        objStar = compute_obj(Xstar, Ystar, problem, (VSPFLOAT) gamma) ;

        VSPINT imbStar = IMBALANCE(Xstar->weight, Ystar->weight);

        if (objStar > objBest || (objStar == objBest && imbStar < imbBest))
        {
            foundBetter = true;
            objBest = objStar;
            imbBest = imbStar;

            if (Xbest != NULL) Destroy_State(Xbest);
            if (Ybest != NULL) Destroy_State(Ybest);
            Xbest = Xstar;
            Ybest = Ystar;
            Xstar = NULL;
            Ystar = NULL;
        }
        else
        {
            Destroy_State(Xstar);
            Destroy_State(Ystar);
        }
    }

    if (foundBetter)
    {
        Swap_Mem(X, Xbest, sizeof(State));
        Swap_Mem(Y, Ybest, sizeof(State));
        Destroy_State(Xbest);
        Destroy_State(Ybest);
    }
}
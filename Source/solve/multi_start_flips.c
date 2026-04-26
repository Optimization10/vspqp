
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

    VSPINT objBest ;
    Compute_ObjValSimple
    (
        X->point->int_part, X->Az->int_part, X->weight,
        Y->point->int_part, Y->Az->int_part, Y->weight,
        problem, &objBest
    ) ;
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

        Greedy_Flips(Xstar, Ystar, problem, numPassesInFlips, true, false);

        VSPINT objStar;
        Compute_ObjValSimple
        (
            Xstar->point->int_part, Xstar->Az->int_part, Xstar->weight,
            Ystar->point->int_part, Ystar->Az->int_part, Ystar->weight,
            problem, &objStar
        );
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
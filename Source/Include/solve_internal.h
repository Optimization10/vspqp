#ifndef SOLVE_INTERNAL_H
#define SOLVE_INTERNAL_H

#include "core.h"

/*----------------------------------------------------------------------------*/
/*  Declarations                                                              */
/*----------------------------------------------------------------------------*/
void set_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
) ;

void Optimize
(
    State                     *X,
    State                     *Y,
    Problem             *problem,
    Options             *options
) ;

void Multi_Start_Flips
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       numFlipRuns,
    VSPINT       numPassesInFlips
) ;


void Best_Of_MC_And_Greedy_Flips
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       numFlipsRuns
) ;


typedef void (*Strategy)(State*, State*, Problem*, Options*) ;

bool Try_And_Keep_If_Better
(
    State *X, State *Y, Problem *problem,
    Strategy strategy, Options *options, VSPFLOAT *ratio 
);

void strategy_mc_and_flips(State*, State*, Problem*, Options*);
void strategy_flips(State*, State*, Problem*, Options*);
void strategy_fm(State*, State*, Problem*, Options*);
void strategy_centroid_restart(State*, State*, Problem*, Options*);
void strategy_bisection_restart
(
    State *X, State *Y, Problem *problem, Options *options
) ;
void strategy_random_restart
(
    State *X, State *Y, Problem *problem, Options *options
) ;
#endif
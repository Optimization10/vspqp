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


void Part_To_States
(
    int        *part,
    State      *X,
    State      *Y,
    Problem    *problem,
    bool        empty_sep
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

void strategy_explore_bisection
(
    State *X, State *Y, Problem *problem, Options *options
) ;

void Rebalance_Singletons(State *X, State *Y, Problem *P) ;

int Cover_Refine(State *X, State *Y, Problem *P, VSPINT max_bip) ;

/*----------------------------------------------------------------------------*/
/*  Multilevel solve pieces shared by VSP_Solve and VSP_Solve_2 (solve.c)     */
/*----------------------------------------------------------------------------*/
// Improves a solution at one level of the hierarchy
typedef void (*Level_Solver)
(
    State    *X,
    State    *Y,
    Problem  *problem,
    Options  *options
) ;

// The full level solve: Optimize, the optional second start, and
// Rebalance_Singletons
void Solve_Level(State *X, State *Y, Problem *problem, Options *options) ;

// Solves the problem hierarchy below probFinest and returns the solution at
// the finest level in (X_0, Y_0)
typedef void (*Hierarchy_Solver)
(
    Problem  *probFinest,
    Options  *options,
    State   **X_0,
    State   **Y_0
) ;

solve_result *VSP_Solve_With
(
    VSPINT             n_input,
    VSPLONG           *Ap_input,
    VSPINT            *Ai_input,
    Options           *options,
    char              *name,
    Hierarchy_Solver   solve_hierarchy
) ;

bool Solve_Coarsest
(
    Problem      *P,
    Options      *options,
    VSPINT       *sol_sings,
    VSPINT        n_sings,
    Level_Solver  solve_level,
    State       **X,
    State       **Y
) ;

void Solve_Up
(
    Problem      *P,
    Problem      *Top,
    bool          have_global_optimum,
    Options      *options,
    Level_Solver  solve_level,
    State       **X,
    State       **Y
) ;

void Destroy_All_Problems(Problem *problemFinest) ;
#endif
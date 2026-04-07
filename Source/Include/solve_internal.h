#ifndef SOLVE_INTERNAL_H
#define SOLVE_INTERNAL_H

#include "core.h"

/*----------------------------------------------------------------------------*/
/*  Declarations                                                              */
/*----------------------------------------------------------------------------*/

void ResolveLevelDecideAndDo
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options
) ;

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

void Copy_Data
(
    VSPINT *xfrom,
    VSPINT *Axfrom,
    VSPINT *wgtxfrom,
    VSPINT *yfrom,
    VSPINT *Ayfrom,
    VSPINT *wgtyfrom,
    VSPINT *xto,
    VSPINT *Axto,
    VSPINT *wgtxto,
    VSPINT *yto,
    VSPINT *Ayto,
    VSPINT *wgtyto,
    VSPINT n
) ;

void Best_Of_MC_And_Greedy_Flips
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options,
    VSPINT  numFlipsRuns
) ;


typedef void (*Strategy)(State*, State*, Problem*, Options*);

bool Try_And_Keep_If_Better
(
    State *X, State *Y, Problem *problem,
    Strategy strategy, Options *options
);

void strategy_mc_and_flips(State*, State*, Problem*, Options*);
void strategy_flips(State*, State*, Problem*, Options*);
void strategy_fm(State*, State*, Problem*, Options*);
void strategy_centroid_restart(State*, State*, Problem*, Options*);

#endif
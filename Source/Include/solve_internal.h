#ifndef SOLVE_INTERNAL_H
#define SOLVE_INTERNAL_H

#include "shared_api.h"
#include "coarsen_api.h"

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

void Optimize
(
    State                     *X,
    State                     *Y,
    Problem             *problem,
    Options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
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

#endif
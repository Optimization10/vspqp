#ifndef FLIPS_PUBLIC_V5_H
#define FLIPS_PUBLIC_V5_H
#include "shared_api.h"


void Flips_v5
(
    VSPINT                     *x,
    VSPINT                    *Ax,
    VSPINT                  *wgtx,
    VSPINT                     *y,
    VSPINT                    *Ay,
    VSPINT                  *wgty,
    Problem          *problem,
    VSPINT           maxNumPasses
) ;

void TryFlipsManyStartingPoints
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
#endif
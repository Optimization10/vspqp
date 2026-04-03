#ifndef FLIPS_PUBLIC_V5_H
#define FLIPS_PUBLIC_V5_H
#include "shared_api.h"

// WARNING: If this function is called in debugging mode, then, since it uses
// the random number generator, it will make the results different from when
// running in performance mode
void Flips_v5_UnitTest
(
    Problem *problem
) ;

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
    VSP_options               *options,
    VSPINT                 numFlipRuns,
    VSPINT            numPassesInFlips
) ;
#endif
#ifndef FLIPS_PUBLIC_V7_H
#define FLIPS_PUBLIC_V7_H
#include "shared_api.h"

// WARNING: If this function is called in debugging mode, then, since it uses
// the random number generator, it will make the results different from when
// running in performance mode
void Flips_v7_UnitTest
(
    Problem *problem
) ;

void Flips_v7
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

#endif
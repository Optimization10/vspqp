#ifndef FLIPS_PUBLIC_V4_H
#define FLIPS_PUBLIC_V4_H
#include "shared_api.h"


void Flips_v4
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
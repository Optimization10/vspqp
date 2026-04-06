#ifndef FLIPS_PUBLIC_V4_H
#define FLIPS_PUBLIC_V4_H
#include "shared_api.h"
#include <stdbool.h>

void Flips_Dep
(
    State               *X,
    State               *Y,
    Problem       *problem,
    VSPINT        maxNumPasses
) ;



void Flips_v999
(
    State *X,
    State *Y,
    Problem              *problem,
    VSPINT           maxNumPasses,
    bool              deep_search
) ;

#endif
#ifndef FLIPS_PUBLIC_V4_H
#define FLIPS_PUBLIC_V4_H
#include "shared_api.h"
#include <stdbool.h>

void Flips
(
    State               *X,
    State               *Y,
    Problem       *problem,
    VSPINT        maxNumPasses
) ;



void Greedy_Flips
(
    State *X,
    State *Y,
    Problem              *problem,
    VSPINT           maxNumPasses,
    bool              deep_search
) ;

#endif
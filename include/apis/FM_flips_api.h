#ifndef FM_FLIPS_API_H
#define FM_FLIPS_API_H

#include "shared_api.h"

void FM_Flips_API_v2
(
    VSP_state         *Z,
    VSP_state         *W,
    VSP_problem *problem,
    VSP_options *options
) ;

void FM_Flips_v7
(
    VSP_state              *X,
    VSP_state              *Y,
    VSP_problem      *problem,
    VSPINT       maxNumPasses
) ;

void FM_Flips_v8
(
    VSP_state              *X,
    VSP_state              *Y,
    VSP_problem      *problem,
    VSPINT       maxNumPasses
) ;

void FM_Flips_v9
(
    VSP_state               *X,
    VSP_state               *Y,
    VSP_problem       *problem,
    VSPINT        maxNumPasses
) ;

#endif
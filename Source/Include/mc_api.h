#ifndef MC_API_H
#define MC_API_H

#include "shared_api.h"

typedef enum
{
    NO_REDUCE,
    USE_MINREDUCE,
    REDUCE_SIMPLY
} GR_Reduce_Type ;


VSPFLOAT Reduce_Gamma
(
    VSPFLOAT       gammaIn,
    GR_Reduce_Type redType,
    State           *Z,
    State           *W,
    Options   *options,
    Problem   *problem
) ;

/* This file contains functions involving mc which are to be exposed to the
   outside */
void VSP_mountain_climbing
(
    State               *X,
    State               *Y,
    Problem       *problem,
    Options       *options,
    VSPFLOAT             gamma,
    bool          is_gamma_org
) ;


void VSP_mountain_climbing_v2
(
    State                  *X,
    State                  *Y,
    Problem          *problem,
    Options          *options,
    VSPFLOAT                gamma,
    VSPINT                 numLPs /* >= 0. If 0, we launch directly into ortho*/
) ;

VSPFLOAT VSP_LP_solver_4
(
    State *Z, /* Solution */
    VSP_Az *Aw,
    VSPINT *c,
    VSPINT u,
    VSPFLOAT gamma,
    Problem *problem,
    Options *options,
    char             var  /* Variable we are optimizing over */
) ;

#endif
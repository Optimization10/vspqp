#ifndef MC_API_H
#define MC_API_H



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
void mountain_climbing
(
    State               *X,
    State               *Y,
    Problem       *problem,
    Options       *options,
    VSPFLOAT             gamma,
    bool          is_gamma_org
) ;


void mountain_climbing_v2
(
    State                  *X,
    State                  *Y,
    Problem          *problem,
    Options          *options,
    VSPFLOAT                gamma,
    VSPINT                 numLPs /* >= 0. If 0, we launch directly into ortho*/
) ;

VSPFLOAT LP_solver_4
(
    State *Z, /* Solution */
    Az *Aw,
    VSPINT *c,
    VSPINT u,
    VSPFLOAT gamma,
    Problem *problem,
    Options *options,
    char             var  /* Variable we are optimizing over */
) ;

#endif
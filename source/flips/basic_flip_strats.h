#ifndef BASIC_FLIP_STRATS_H
#define BASIC_FLIP_STRATS_H

#include "shared_api.h"
#include "flips_context.h"
#include "ds_api.h"


void BasicFlipStrats_FeasUpflipsUntilTargetWgt
(
    FlipsCtx                  *FC,
    VSPINT               whichVar, // Can NOT be both variables
    bool              withLocking,
    Flip_Selector_Mode       mode,
    VSPINT              targetWgt
) ;

void BasicFlipStrats_FeasBeneficialFlips
(
    FlipsCtx                  *FC,
    VSPINT              whichVars,
    bool              withLocking,
    VSPFLOAT                gamma
) ;

#endif
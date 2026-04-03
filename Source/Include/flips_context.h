#ifndef FLIPS_CONTEXT_H
#define FLIPS_CONTEXT_H

#include "shared_api.h"

typedef struct FlipsCtx FlipsCtx ;

typedef struct
{
    VSPINT    ind ;
    VSPINT    var ;
    VSPFLOAT gain ;
} Flip ;

typedef enum
{
    FEAS_FLIP,
    FEAS_BEN_FLIP

} Flip_Selector_Mode ;

typedef enum
{
    VAR_0      = 0,
    VAR_1      = 1,
    BOTH_VARS  = 2
} Flip_Selector_Vars ;

FlipsCtx *FlipsCtx_Init
(
    VSP_state   *states[2],
    bool         randomize,
    Problem   *problem
) ;

void FlipsCtx_Destroy(FlipsCtx *FC) ;

// whichVars = 0 --- var 0
// whichVars = 1 --- var 1
// whichVars = 2 --- both vars

// mode = 0 --- select the best feasible, even with negative gain 
// mode = 1 --- select the best beneficial flip. Return EMPTY if no
//              feas ben flip
Flip FlipsCtx_Select_Flip
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars,
    Flip_Selector_Mode           mode,
    VSPFLOAT                    gamma
) ;

void FlipsCtx_Reset
(
    FlipsCtx *FC
) ;

// Let z be the variable of the flip, and w the other variable
// ASSUMED: The flip is happening with an Awi or Azi with optimal T
void FlipsCtx_Apply_Flip
(
    FlipsCtx      *FC,
    Flip        *flip,
    bool       lockIt  // lock the flip we do
) ;

bool FlipsCtx_Is_Nonempty_Flip(Flip *flip) ;

// If this is an empty flip, return infeasible
// If u is empty, then problem->u is used
bool FlipsCtx_Is_Feasible_Flip(FlipsCtx *FC, Flip *flip, VSPINT u) ;


void FlipsCtx_Check
(
    FlipsCtx *FC
) ;

VSPINT FlipsCtx_Get_Wgt
(
    FlipsCtx  *FC,
    VSPINT    var
) ;

void FlipsCtx_Set_Wgt
(
    FlipsCtx  *FC,
    VSPINT    var,
    VSPINT    wgt
) ;

void FlipsCtx_Initialize_Lists_And_Ts(FlipsCtx *FC) ;

VSPINT FlipsCtx_GetFlipDir
(
    FlipsCtx    *FC,
    Flip      *flip
) ;


Flip FlipsCtx_Select_Flip_FMv9
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars
) ;

// Completely ignore everything except vars, Avars, and wgts
void FlipsCtx_Apply_Flip_Basic
(
    FlipsCtx      *FC,
    Flip        *flip
) ;

#endif
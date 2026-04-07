#ifndef UTILS_H
#define UTILS_H

#include "core.h"

State *Clone_State(State *src) ;

void Round_Down_Fractions
(
    State           *Z,
    State           *W,
    Problem   *problem
) ;

void Swap_Mem(void *a, void *b, size_t size) ;

void Az_pm_col_of_A
(
    VSPINT             *Az, // vector to which we add/subtract
    VSPINT               i, // column to add/subtract
    VSPINT          operation, // 1 if add, -1 if subtract
    Problem   *problem
) ;

void Compute_Az
(
    Point              *z,
    Az                *Az,
    VSPINT                *Ai,
    VSPLONG               *Ap,
    VSPINT                  n
) ;

#endif

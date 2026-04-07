#ifndef MC_INTERNAL_H
#define MC_INTERNAL_H

#include "shared_api.h"
#include "mc_api.h" // IWYU pragma: keep


/******************************************************************************
    TYPES
*******************************************************************************/


typedef struct value_index_struct
{
    VSPFLOAT float_value;
    VSPINT int_value;
    VSPINT index;
} value_index ;




/******************************************************************************
    FUNCTIONS
*******************************************************************************/

VSPFLOAT LP_solver_0
(
    VSPFLOAT          *x, /* solution                                         */
    VSPFLOAT          *d, /* objective                                        */
    VSPINT            *c, /* constraint                                       */
    VSPINT             l, /* Lower bound. Must be 0 or 1                      */
    VSPINT             u, /* RHS of constraint. Must be at least 1            */
    VSPINT             n, /* size of problem                                  */
    Options *options  /* Choose randomness or not                         */
) ;

#endif
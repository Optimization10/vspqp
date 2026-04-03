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

VSPFLOAT LP_solver_0 (
    VSPFLOAT          *x, /* solution                                         */
    VSPFLOAT          *d, /* objective                                        */
    VSPINT            *c, /* constraint                                       */
    VSPINT             l, /* Lower bound. Must be 0 or 1                      */
    VSPINT             u, /* RHS of constraint. Must be at least 1            */
    VSPINT             n, /* size of problem                                  */
    Options *options  /* Choose randomness or not                         */
) ;

bool Is_LP_Optimal
(
  VSPINT              *z,
  VSPINT               k,
  VSPFLOAT            zk,
  VSPINT           wgt_z,
  Az             *Aw,
  VSPINT              *c,
  VSPINT               u,
  VSPFLOAT         gamma,
  VSPINT      sum_all_ci,
  Options   *options
) ;

bool LP_is_optimal
(
    State         *Z,
    Az           *Aw,
    VSPFLOAT       gamma,
    Problem *problem,
    Options *options
) ;


bool LP_is_optimal_case_1a
(
    Point          *z,
    Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool LP_is_optimal_case_1b
(
    Point          *z,
    Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool LP_is_optimal_case_2a
(
    Point          *z,
    Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool LP_is_optimal_case_2b
(
    Point          *z,
    Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool LP_is_optimal_case_3
(
    Point          *z,
    Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

#endif
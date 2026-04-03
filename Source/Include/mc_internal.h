#ifndef MC_INTERNAL_H
#define MC_INTERNAL_H

#include "shared_api.h"
#include "mc_api.h" // IWYU pragma: keep


/******************************************************************************
    TYPES
*******************************************************************************/


typedef struct VSP_value_index_struct
{
    VSPFLOAT float_value;
    VSPINT int_value;
    VSPINT index;
} VSP_value_index ;




/******************************************************************************
    FUNCTIONS
*******************************************************************************/

VSPFLOAT VSP_LP_solver_0 (
    VSPFLOAT          *x, /* solution                                         */
    VSPFLOAT          *d, /* objective                                        */
    VSPINT            *c, /* constraint                                       */
    VSPINT             l, /* Lower bound. Must be 0 or 1                      */
    VSPINT             u, /* RHS of constraint. Must be at least 1            */
    VSPINT             n, /* size of problem                                  */
    VSP_options *options  /* Choose randomness or not                         */
) ;

bool Is_LP_Optimal
(
  VSPINT              *z,
  VSPINT               k,
  VSPFLOAT            zk,
  VSPINT           wgt_z,
  VSP_Az             *Aw,
  VSPINT              *c,
  VSPINT               u,
  VSPFLOAT         gamma,
  VSPINT      sum_all_ci,
  VSP_options   *options
) ;

bool VSP_LP_is_optimal
(
    VSP_state         *Z,
    VSP_Az           *Aw,
    VSPFLOAT       gamma,
    Problem *problem,
    VSP_options *options
) ;


bool VSP_LP_is_optimal_case_1a
(
    Point          *z,
    VSP_Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool VSP_LP_is_optimal_case_1b
(
    Point          *z,
    VSP_Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool VSP_LP_is_optimal_case_2a
(
    Point          *z,
    VSP_Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool VSP_LP_is_optimal_case_2b
(
    Point          *z,
    VSP_Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

bool VSP_LP_is_optimal_case_3
(
    Point          *z,
    VSP_Az            *Aw,
    VSPFLOAT        gamma,
    VSPINT             *c
) ;

#endif
#include "core.h"
#include "mc_internal.h"
#include "utils.h"
#include <stdio.h>


VSPFLOAT minreduce_one_var
(
    State         *Z,
    Az           *Aw,
    Problem *problem
) ;


/* Return the smallest gamma_t in (0, gamma) such that (x, y) is
   a stationary point                                                         */
VSPFLOAT minreduce
(
    State         *X,
    State         *Y,
    Options *options,
    Problem *problem
)
{

    static long long unsigned count = 0 ;
    count++ ;

    VSPFLOAT gamma_tx = minreduce_one_var(X, Y->Az, problem) ;
    VSPFLOAT gamma_ty = minreduce_one_var(Y, X->Az, problem) ;
    VSPFLOAT gamma_max = (gamma_tx > gamma_ty) ? gamma_tx : gamma_ty ;

    #ifdef DEBUG
        check_states(X, Y, problem) ;
    #endif

    return gamma_max ;
}

VSPFLOAT minreduce_one_var
(
    State         *Z,
    Az           *Aw,
    Problem *problem
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    if(Z->weight == 0 || Z->weight == problem->u) {
        return EMPTY ;
    }

    VSPINT    i ;
    VSPINT    n              = problem->graph->n ;
    VSPINT    one            = Aw->one ;
    VSPINT   *Aw_int         = Aw->int_part ;
    VSPINT   *col_A          = Aw->col_frac_z ;
    VSPINT   *z_int          = Z->point->int_part ;
    VSPINT   *c              = problem->graph->c ;
    VSPFLOAT  val_frac_w     = Aw->val_frac_z ;
    VSPFLOAT  Awi ;
    VSPFLOAT  gamma_t        = EMPTY ;
    VSPFLOAT  gamma_second   = EMPTY ;
    bool      w_has_frac_ind = (Aw->ind_frac_z != EMPTY) ;

    for(i = 0 ; i < n ; i++)
    {
        Awi = Aw_int[i] + w_has_frac_ind * (col_A[i] == one) * val_frac_w ;
        if (z_int[i] == 0 && Awi > 0)
        {
            if(c[i] / Awi > gamma_t)
            {
                gamma_second = gamma_t ;
                gamma_t = c[i] / Awi ;
            }
        }
    }


    if(gamma_second != EMPTY) return gamma_second ;
    else                      return gamma_t ;
}
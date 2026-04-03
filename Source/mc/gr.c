#include "shared_api.h"
#include "gr_api.h"
#include "mc_api.h"



VSPFLOAT VSP_minreduce_one_var
(
    State         *Z,
    VSP_Az           *Aw,
    Problem *problem
) ;


/* Return the smallest gamma_t in (0, gamma) such that (x, y) is
   a stationary point                                                         */
VSPFLOAT VSP_minreduce
(
    State         *X,
    State         *Y,
    VSP_options *options,
    Problem *problem
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPFLOAT gamma_tx = VSP_minreduce_one_var(X, Y->Az, problem) ;
    VSPFLOAT gamma_ty = VSP_minreduce_one_var(Y, X->Az, problem) ;

    VSPFLOAT gamma_max = (gamma_tx > gamma_ty) ? gamma_tx : gamma_ty ;


    if( gamma_max >= problem->Graph->max_vtx_cost )
    {
        VSP_mountain_climbing(X, Y, problem, options, problem->Graph->max_vtx_cost, true) ;
    }
    gamma_tx = VSP_minreduce_one_var(X, Y->Az, problem) ;
    gamma_ty = VSP_minreduce_one_var(Y, X->Az, problem) ;

    gamma_max = (gamma_tx > gamma_ty) ? gamma_tx : gamma_ty ;

    if(gamma_max > problem->Graph->max_vtx_cost)
    {
        // printf("Reduced gamma greater. Problem: %s\n", problem->name) ;
    }
    else if(gamma_max < 0)
    {
        // printf("Negative\n") ;
    }

    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
    #endif

    return gamma_max ;
}

VSPFLOAT VSP_minreduce_one_var
(
    State         *Z,
    VSP_Az           *Aw,
    Problem *problem
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    if(Z->weight == 0 || Z->weight == problem->u) {
        return EMPTY ;
    }

    VSPINT    i ;
    VSPINT    n              = problem->Graph->n ;
    VSPINT    one            = Aw->one ;
    VSPINT   *Aw_int         = Aw->int_part ;
    VSPINT   *col_A          = Aw->col_frac_z ;
    VSPINT   *z_int          = Z->point->int_part ;
    VSPINT   *c              = problem->Graph->c ;
    VSPFLOAT  val_frac_w     = Aw->val_frac_z ;
    VSPFLOAT  Awi ;
    VSPFLOAT  gamma_t        = EMPTY ;
    VSPFLOAT  gamma_second   = EMPTY ;
    bool      w_has_frac_ind = (Aw->ind_frac_z != EMPTY) ;

    for(i = 0 ; i < n ; i++)
    {
        Awi = Aw_int[i] + w_has_frac_ind * (col_A[i] == one) * val_frac_w ;
        if (z_int[i] == 0 && Awi > 0) {
            if(c[i] / Awi > gamma_t)
            {
                gamma_second = gamma_t ;
                gamma_t = c[i] / Awi ;
                // if(gamma_t == 4 && problem->Graph->max_vtx_cost == 2) printf("\ni: %lld", (long long) i) ;
            }
        }
        #ifdef VSP_DEBUG
            if(z_int[i] == 0 && Awi == 0)
            {
                // Error("Entered gamma reduction when not neccesary", DBG_INFO) ;
            }
        #endif
    }

    #ifdef VSP_DEBUG
        if(gamma_t > problem->Graph->max_vtx_cost) {
            Error("Error in gamma reduction", DBG_INFO) ;
        }
    #endif

    if(gamma_second != EMPTY) return gamma_second ;
    return gamma_t ;
}


static void i_am_dead_code(void)
{
    int x = 1;
}

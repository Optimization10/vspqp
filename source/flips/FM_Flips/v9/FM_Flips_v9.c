#include "FM_Flips_v9.h"
#include "shared_api.h"
#include "math.h"
#include "score.h"
#include "states_manager.h"
#include "basic_flip_strats.h"
#include "flips_mem.h"


void FM_Flips_v9
(
    VSP_state               *X,
    VSP_state               *Y,
    VSP_problem       *problem,
    VSPINT        maxNumPasses
)
{
    #ifdef VSP_DEBUG
        VSP_check_states(X, Y, problem) ;
        VSPFLOAT objIn = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
    #endif
    SET_FUNCTION_COUNTER ;
    bool lock ;
    VSPINT dir,
    n = problem->Graph->n ;
    FlipsCtx *FC = FlipsCtx_Init((VSP_state *[2]){X, Y}, true, problem);
    FlipsMemory *FM = FlipsMem_Construct(12 * n, FC) ;

    for(VSPINT pass = 1 ; pass <= maxNumPasses ; pass++)
    {
        #ifdef VSP_DEBUG
            VSP_check_states(X, Y, problem) ;
            VSPFLOAT objIn = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        #endif

        // Initialize structs
        if(pass > 1)
        {
            FlipsMem_Reset(FM) ;
            FlipsCtx_Reset(FC) ;
        }

        // Do flips
        while( 1 )
        {
            Flip flip = FlipsCtx_Select_Flip_FMv9(FC, 2) ;
            if( !FlipsCtx_Is_Nonempty_Flip(&flip) ) break ;
            dir = FlipsCtx_GetFlipDir(FC, &flip) ;
            lock = (dir == 1) ;
            FlipsCtx_Apply_Flip(FC, &flip, lock) ;
            FlipsMem_Update(FM, &flip) ;
        }

        // Roll back to best point
        FlipsMem_Rollback(FM) ;
        X->weight = FlipsCtx_Get_Wgt(FC, 0) ;
        Y->weight = FlipsCtx_Get_Wgt(FC, 1) ;

        #ifdef VSP_DEBUG
            VSP_check_states(X, Y, problem) ;
            VSPFLOAT objFin = VSP_compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
            if(!VSP_equal(objFin, objIn + FlipsMem_GetBestAccumGain(FM)) )
            {
                ERR("Obj did not increase as expected") ;
            }
        #endif
    }



    FlipsMem_Destroy(FM) ;
    FlipsCtx_Destroy(FC) ;


    return ;
}
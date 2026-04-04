#include "core.h"
#include "flips_private_v5.h"
#include "tmp.h"

VSPINT g_num_vert_curr_lvl_v5 ;

void Flips_v5
(
    VSPINT                     *x,
    VSPINT                    *Ax,
    VSPINT                  *wgtx,
    VSPINT                     *y,
    VSPINT                    *Ay,
    VSPINT                  *wgty,
    Problem          *problem,
    VSPINT           maxNumPasses
)
{
    SET_FUNCTION_COUNTER ;
    g_num_vert_curr_lvl_v5 = problem->graph->n ;
    #ifdef DEBUG
        Check_Inputs_v5(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
    #endif


    VSPINT n     = problem->graph->n ;
    VSPINT *c    = problem->graph->c ;
    VSPINT u     = problem->u ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    /*------------------------------------------------------------------------*/
    // Initialize all structs
    /*------------------------------------------------------------------------*/
    VSPINT *maxGroupSizes = malloc((gamma + 1) * sizeof(*maxGroupSizes)) ;
    Compute_Max_Group_Sizes_v5(maxGroupSizes, c, gamma, n) ;

    // x structs
    IdxGroupsByC *x0Ay0_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 0, 0, x, Ay, 'x', MAX_HEAP_FLAG) ;
    Direct_Access_List *x0Ay1_list = Construct_Direct_Access_List(n) ;

    IdxGroupsByC *x1Ay1_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 1, 1, x, Ay, 'x', MIN_HEAP_FLAG) ;
    IdxGroupsByC *x1Ay2_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 1, 2, x, Ay, 'x', MIN_HEAP_FLAG) ;


    // y structs
    IdxGroupsByC *y0Ax0_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 0, 0, y, Ax, 'y', MAX_HEAP_FLAG) ;
    Direct_Access_List *y0Ax1_list = Construct_Direct_Access_List(n) ;

    IdxGroupsByC *y1Ax1_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 1, 1, y, Ax, 'y', MIN_HEAP_FLAG) ;
    IdxGroupsByC *y1Ax2_groups = IdxGroupsByC_Init(gamma, maxGroupSizes, c, n, 1, 2, y, Ax, 'y', MIN_HEAP_FLAG) ;

    Initialize_All_Structs_v5
    (
        x, Ax, *wgtx, y, Ay, *wgty,
        x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
        y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
        problem
    ) ;

    VSPINT numFeasibleStrictlyBenFlips, numFeasibleBenFlips ;

    FeasibleBeneficialFlips_v5
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        problem, &numFeasibleBenFlips
    ) ;

    
    Neutral_Flips_Strategy_1_v5
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        maxNumPasses, problem
    ) ;

    Balance_Flips
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        problem
    ) ;
    
    FeasibleBeneficialFlips_v5
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        problem, &numFeasibleBenFlips
    ) ;


    // Free memory
    free(maxGroupSizes) ;
    IdxGroupsByC_Destroy(x0Ay0_groups) ;
    IdxGroupsByC_Destroy(x1Ay1_groups) ;
    IdxGroupsByC_Destroy(x1Ay2_groups) ;

    IdxGroupsByC_Destroy(y0Ax0_groups) ;
    IdxGroupsByC_Destroy(y1Ax1_groups) ;
    IdxGroupsByC_Destroy(y1Ax2_groups) ;

    Destruct_Direct_Access_List(x0Ay1_list) ;
    Destruct_Direct_Access_List(y0Ax1_list) ;
}

void StrictlyBeneficialFlips_v5
(
    VSPINT                            *x,
    VSPINT                            *y,
    VSPINT                           *Ax,
    VSPINT                           *Ay,
    VSPINT                         *wgtx,
    VSPINT                         *wgty,

    IdxGroupsByC           *x0Ay0_groups,
    IdxGroupsByC           *x1Ay1_groups,
    IdxGroupsByC           *x1Ay2_groups,
    Direct_Access_List       *x0Ay1_list,

    IdxGroupsByC           *y0Ax0_groups,
    IdxGroupsByC           *y1Ax1_groups,
    IdxGroupsByC           *y1Ax2_groups,
    Direct_Access_List       *y0Ax1_list,

    Problem                 *problem,
    VSPINT              *retNumFlipsDone
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        VSPINT objDbgInit ;
        Compute_ObjValDbg(x, Ax, y, Ay, problem, &objDbgInit) ;
    #endif


    /*------------------------------------------------------------------------*/
    // Do flips
    /*------------------------------------------------------------------------*/
    char flipVar ;
    VSPINT indFlip ;
    VSPINT numFlips = 0 ;
    while( 1 )
    {
        // Choose the best flip
        SelectRandomOptimalStrictlyBeneficialFlip_v5
        (
            x, y, Ax, Ay, *wgtx, *wgty,
            x0Ay0_groups,
            x1Ay1_groups,
            x1Ay2_groups,
            x0Ay1_list,
        
            y0Ax0_groups,
            y1Ax1_groups,
            y1Ax2_groups,
            y0Ax1_list,
        
            problem,
            &flipVar, // '\0' if no flip
            &indFlip  // EMPTY if no flip
        ) ;

        #ifdef LOG_FLIPS_V5
            if( LOG_CONDITION )
            {
                if(indFlip != EMPTY)
                {
                    char filename[150] ;
                    snprintf(filename, sizeof(filename), "Flips_v5/problem_%s_Level_%lld.csv", problem->name, (long long) problem->level) ;
                    Log(filename, "%c, %lld", flipVar, (long long) indFlip) ;
                }
            }
        #endif

        // Do the best flip
        if(flipVar == 'x')
        {
            Do_Flip
            (
                x, Ax, wgtx, y, Ay, wgty,
                x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
                y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
                indFlip, problem
            ) ;
        }
        else if(flipVar == 'y')
        {
            Do_Flip
            (
                y, Ay, wgty, x, Ax, wgtx,
                y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
                x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
                indFlip, problem
            ) ;
        }
        else
        {
            break ;
        }
        numFlips++ ;
    }

    #ifdef DEBUG
        Check_Outputs_v5(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
        VSPINT objDbgFinal ;
        Compute_ObjValDbg(x, Ax, y, Ay, problem, &objDbgFinal) ;
        if( !(objDbgFinal >= objDbgInit) )
        {
            ERR("Objective was damaged in strictly beneficial feasible flips") ;
        }
    #endif

    // Prepare return
    *retNumFlipsDone = numFlips ;
    return ;
}
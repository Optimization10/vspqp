#include "core.h"
#include "greedy_flips_internal.h"
#include <stdbool.h>

void FeasibleBeneficialFlips
(
    VSPINT                            *x,
    VSPINT                            *y,
    VSPINT                           *Ax,
    VSPINT                           *Ay,
    VSPINT                         *wgtx,
    VSPINT                         *wgty,
    bool                     deep_search,
    bool                          strict,

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
        if(!strict) // Not sure why I'm guarding, just copying what I did before
        {
            Check_All_Data
            (
                x, Ax, *wgtx, y, Ay, *wgty,
                x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
                y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
                problem
            ) ;
        }

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
        SelectRandomOptimalBeneficialFlip
        (
            x, y, Ax, Ay, *wgtx, *wgty, deep_search, strict,
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
        if(deep_search)
        {
            Check_Outputs_v5(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
        }
        else
        {
            Check_Outputs(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
        }
        
        VSPINT objDbgFinal ;
        Compute_ObjValDbg(x, Ax, y, Ay, problem, &objDbgFinal) ;
        if( !(objDbgFinal >= objDbgInit) )
        {
            ERR("Objective was damaged in strictly beneficial feasible flips") ;
        }

        // Check that there are no neutral downflips available
        if(!strict)
        {
            for(VSPINT i = 0 ; i < problem->graph->n ; i++)
            {
                if(x[i] && Ay[i] && *wgtx > problem->graph->c[i])
                {
                    ERR("Downflip still available") ;
                }

                if(y[i] && Ax[i] && *wgty > problem->graph->c[i])
                {
                    ERR("Downflip still available") ;
                }
            }
        }
    #endif

    // Prepare return
    *retNumFlipsDone = numFlips ;
    return ;
}

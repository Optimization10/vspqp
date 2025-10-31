#include "flips_private_v4.h"

void FeasibleBeneficialFlips
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

    VSP_problem                 *problem,
    VSPINT              *retNumFlipsDone
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        Check_All_Data
        (
            x, Ax, *wgtx, y, Ay, *wgty,
            x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
            y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
            problem
        ) ;

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

        #ifdef LOG_FLIPS_V4
            // Log the flip
            if( LOG_CONDITION )
            {
                if(indFlip != EMPTY)
                {
                    char filename[150] ;
                    snprintf(filename, sizeof(filename), "Flips_v4/problem_%s_Level_%lld.csv", problem->name, (long long) problem->level) ;
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

    #ifdef VSP_DEBUG
        Check_Outputs(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
        VSPINT objDbgFinal ;
        Compute_ObjValDbg(x, Ax, y, Ay, problem, &objDbgFinal) ;
        if( !(objDbgFinal >= objDbgInit) )
        {
            ERR("Objective was damaged in strictly beneficial feasible flips") ;
        }

        // Check that there are no neutral downflips available
        for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
        {
            if(x[i] && Ay[i] && *wgtx > problem->Graph->c[i])
            {
                ERR("Downflip still available") ;
            }

            if(y[i] && Ax[i] && *wgty > problem->Graph->c[i])
            {
                ERR("Downflip still available") ;
            }
        }
    #endif

    // Prepare return
    *retNumFlipsDone = numFlips ;
    return ;
}


void StrictlyBeneficialFlips
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

    VSP_problem                 *problem,
    VSPINT              *retNumFlipsDone
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
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
        SelectRandomOptimalStrictlyBeneficialFlip
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

        #ifdef LOG_FLIPS_V4
            if( LOG_CONDITION )
            {
                if(indFlip != EMPTY)
                {
                    char filename[150] ;
                    snprintf(filename, sizeof(filename), "Flips_v4/problem_%s_Level_%lld.csv", problem->name, (long long) problem->level) ;
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

    #ifdef VSP_DEBUG
        Check_Outputs(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
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
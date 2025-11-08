#include "flips_private_v4.h"


void Flips_v4
(
    VSPINT                     *x,
    VSPINT                    *Ax,
    VSPINT                  *wgtx,
    VSPINT                     *y,
    VSPINT                    *Ay,
    VSPINT                  *wgty,
    VSP_problem          *problem,
    VSPINT           maxNumPasses
)
{
    static long long unsigned count = 0 ;
    count++ ;

    #ifdef VSP_DEBUG
        Check_Inputs(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
    #endif

    /*------------------------------------------------------------------------*/
    static VSPINT nums = 0 ;
    nums++; 
    if(nums == 1)
    {
        printf("\n\nWARNING: Doing a debugging check in performance mode") ;
        fflush(stdout) ;
    }

    for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
    {
        if(x[i] && Ay[i] >= 3)
        {
            printf("\nError: xi is 1 and Ayi is >= 3 in flip. Problem: %s", problem->name) ;
            fflush(stdout) ;
        }

        if(y[i] && Ax[i] >= 3)
        {
            printf("\nError: yi is 1 and Axi is >= 3 in flips. Problem: %s", problem->name) ;
            fflush(stdout) ;
        }
    }
    /*------------------------------------------------------------------------*/


    VSPINT n     = problem->Graph->n ;
    VSPINT *c    = problem->Graph->c ;
    VSPINT u     = problem->u ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;

    /*------------------------------------------------------------------------*/
    // Initialize all structs
    /*------------------------------------------------------------------------*/
    VSPINT *maxGroupSizes = malloc((gamma + 1) * sizeof(*maxGroupSizes)) ;
    Compute_Max_Group_Sizes(maxGroupSizes, c, gamma, n) ;

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

    Initialize_All_Structs
    (
        x, Ax, *wgtx, y, Ay, *wgty,
        x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
        y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
        problem
    ) ;

    VSPINT numFeasibleStrictlyBenFlips, numFeasibleBenFlips ;

    FeasibleBeneficialFlips
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        problem, &numFeasibleBenFlips
    ) ;

    
    Neutral_Flips_Strategy_1
    (
        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        maxNumPasses, problem
    ) ;
    
    FeasibleBeneficialFlips
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

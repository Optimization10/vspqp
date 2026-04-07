#include "core.h"
#include "data_structures.h"
#include "greedy_flips_internal.h"
#include <stdbool.h>
#include "utils.h"



void Greedy_Flips
(
    State      *X,
    State      *Y,
    Problem    *problem,
    VSPINT      maxNumPasses,
    bool        deep_search
)
{
    Round_Down_Fractions(X, Y, problem) ;

    static long long unsigned count = 0 ;
    count++ ;

    VSPINT *x = X->point->int_part ;
    VSPINT *Ax = X->Az->int_part ;
    VSPINT *wgtx = &(X->weight) ;

    VSPINT *y = Y->point->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT *wgty = &(Y->weight) ;

    VSPINT n     = problem->graph->n ;
    VSPINT *c    = problem->graph->c ;
    VSPINT u     = problem->u ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

#ifdef DEBUG
    Check_Inputs(x, Ax, *wgtx, y, Ay, *wgty, problem) ;
#endif

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
        x, y, Ax, Ay, wgtx, wgty, deep_search, false,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
        problem, &numFeasibleBenFlips
    ) ;


    Neutral_Flips_Strategy_1
    (
        x, y, Ax, Ay, wgtx, wgty, deep_search,
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

    FeasibleBeneficialFlips
    (
        x, y, Ax, Ay, wgtx, wgty, deep_search, false,
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


// Return a random, feasible upflip in z with non-negative, highest gain
// May not catch all upflips 
static void Get_Upflip
(
    VSPINT                      wgt_z,
    IdxGroupsByC        *z0Aw0_groups,
    Direct_Access_List    *z0Aw1_list,
    Problem                    *P,
    VSPINT                  *ind_flip
)
{
    VSPINT u = P->u ;

    // Get the maximum ci corresponding to zi = 0, Awi = 0 ;
    VSPINT c_top = Max_Heap_See_Top_Key_v4(z0Aw0_groups->maxHeap) ;

    // If there exists an i with zi = 0, Awi = 0, then randomly select such
    // an upflip and return it
    if(c_top != EMPTY)
    {
        // If not feasible return. No need to check for neutral upflip as
        // that wouldn't be feasible anyway
        if(!FEASIBLE_UPFLIP(wgt_z, c_top, u))
        {
            *ind_flip = EMPTY ;
            return ;
        }

        // Get a random flip
        *ind_flip = Get_Rand_Array_Entry
        (
            z0Aw0_groups->groups[c_top],
            z0Aw0_groups->sizes [c_top]
        ) ;
        return ;
    }

    // No strictly beneficial upflip. Try to return a feasible neutral upflip
    if(!FEASIBLE_UPFLIP(wgt_z, P->graph->max_vtx_cost, u))
    {
        *ind_flip = EMPTY ;
        return ;
    }

    *ind_flip = Direct_Access_List_Get_Rand_Entry(z0Aw1_list) ;
    return ;
}

// Return a downflip with highest gain (all downflips are feasible)
// Assume: z_i = 1 implies Aw_i <= 2
static void Get_Downflip
(
    IdxGroupsByC        *z1Aw2_groups,
    IdxGroupsByC        *z1Aw1_groups,
    Problem                    *P,
    VSPINT                  *ind_flip
)
{
    VSPINT c_bot = EMPTY ;
    IdxGroupsByC *Chosen = NULL ;

    c_bot = Min_Heap_See_Top_Key_v4(z1Aw2_groups->minHeap) ;
    if(c_bot != EMPTY)
    {
        Chosen = z1Aw2_groups ;
    }
    else
    {
        c_bot = Min_Heap_See_Top_Key_v4(z1Aw1_groups->minHeap) ;
        if(c_bot != EMPTY)
        {
            Chosen = z1Aw1_groups ;
        }
    }

    if(c_bot == EMPTY)
    {
        *ind_flip = EMPTY ;
        return ;
    }

    *ind_flip = Get_Rand_Array_Entry
    (
        Chosen->groups[c_bot],
        Chosen->sizes [c_bot]
    ) ;
    return ;
}

// Assume: Input point is the output of feasible beneficial flips
// Attempt to balance the point as much as possible using neutral flips
// and singletons
// Alternate between
// 1) neutral upflip or strictly beneficial upflip in lighter variable
// 2) beneficial downflip in heavier variable
// until the lighter variable surpasses the heavier
// I think this should be called after the last feasible beneficial flips
void Balance_Flips
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

    Problem                       *P
)
{
    // 1. Let z be the lighter variable and w the heavier one
    VSPINT *z;
    VSPINT *w;
    VSPINT *Az;
    VSPINT *Aw;
    VSPINT *wgtz;
    VSPINT *wgtw;
    
    IdxGroupsByC *z0Aw0_groups;
    IdxGroupsByC *z1Aw1_groups;
    IdxGroupsByC *z1Aw2_groups;
    Direct_Access_List *z0Aw1_list;
    
    IdxGroupsByC *w0Az0_groups;
    IdxGroupsByC *w1Az1_groups;
    IdxGroupsByC *w1Az2_groups;
    Direct_Access_List *w0Az1_list;
    
    Assign_By_Weight
    (
        &z, &w, &Az, &Aw, &wgtz, &wgtw,
        &z0Aw0_groups, &z1Aw1_groups, &z1Aw2_groups, &z0Aw1_list,
        &w0Az0_groups, &w1Az1_groups, &w1Az2_groups, &w0Az1_list,

        x, y, Ax, Ay, wgtx, wgty,
        x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
        y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list
    ) ;


    // Alternate flips to promote balance
    VSPINT ind_flip ;
    while(true)
    {
        bool any_flip_done = false ;
        // Do upflip in lighter variable
        Get_Upflip(*wgtz, z0Aw0_groups, z0Aw1_list, P, &ind_flip) ;
        if(ind_flip != EMPTY)
        {
            // Do upflip in lighter variable
            Upflip
            (
                z, Az, w, Aw, wgtz, wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                ind_flip, P
            ) ;
            any_flip_done = true ;
        }

        // Do downflip in heavier variable
        Get_Downflip(w1Az2_groups, w1Az1_groups, P, &ind_flip) ;
        if(ind_flip != EMPTY)
        {
            Downflip
            (
                w, Aw, z, Az, wgtw, wgtz,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                ind_flip, P
            ) ;
            any_flip_done = true ;
        }

        // Break if the lighter variable has surpassed the heavier one
        // TODO 2: May want to continue going if we're improving the objective
        if(*wgtz >= *wgtw || !any_flip_done) break ;
    }
}
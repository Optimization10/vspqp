#include "core.h"
/*------------------------------------------------------------------------------
This file is for the following strategy of neutral flips:

Input: a feasible, binary point (x, y) where the sets{i : xi = 1, Ayi >= 3}
       and {i : yi = 1, Axi >= 3} are empty. (most likely, the input will be a
       point where there are no feasible, strictly beneficial flips)



********************************************************************************
Let S = {i : y_i = 1 and Ay_i = 2}
Repeat 5 times
    Let z be the lighter variable and w the heavier variable

    Initialize M = 1
    Repeat until M = 0:
        Do neutral Upflip_v5s in z. If S ever becomes nonempty, do the Downflip_v5(s)
        contained in S, and return to neutral Upflip_v5s z. Keep doing this until
        there are no more neutral Upflip_v5s in z and S is empty (see Rmk1 below)

        Do all strictly beneficial flips in (z, w), and let M be the number
        of strictly beneficial flips done
    End Repeat

    Initialize M = 1
    Repeat until M = 0:
        Do all neutral Downflip_v5s in w

        Do all strictly beneficial flips in (z, w), and let M be the number of
        strictly beneficial flips done
    End Repeat
End Repeat
Do all strictly beneficial, and neutral Downflip_v5s. Stop when there are no
strictly beneficial Downflip_v5s available,

Rmk 1: It is theoretically possible that doing the Downflip_v5 would violate the
       lower balance constraint. This is very unlikely, but I need to deal with
       it in the code. The code does not deal with this case at the moment. It
       will just do the flip and violate the lower balance constraint, which
       may cause other bugs down the line

Rmk 2: When we say "do all feasible, strictly beneficial flips", it is only
       the feasible, strictly beneficial flips visible to the data structures.
       It is possible that due to being close to the balance constraint, there
       ARE feasible, strictly beneficial Upflip_v5s, but they do not correspond to
       the root of the max heap in the IdxGroupsByC data structs, so they are
       not done 
********************************************************************************

------------------------------------------------------------------------------*/


#include "tmp.h"
#include "flips_private_v5.h"

void Neutral_Flips_Strategy_1_v5
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT *wgtx,
    VSPINT *wgty,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list, // only used in debugging

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list, // only used in debugging

    VSPINT          maxNumPasses,
    Problem       *problem
)
{
    SET_FUNCTION_COUNTER ;
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
    

    VSPINT numFeasibleStrictlyBenFlips = 0 ;
    for(VSPINT numPasses = 0 ; numPasses < maxNumPasses ; numPasses++ )
    {
        /*--------------------------------------------------------------------*/
        // Determine the lighter and heavier variable. Let z be the lighter
        // one and w the heavier one
        /*--------------------------------------------------------------------*/
        // Let z be the lighter variable and w the heavier one
        Assign_By_Weight
        (
            &z, &w, &Az, &Aw, &wgtz, &wgtw,
            &z0Aw0_groups, &z1Aw1_groups, &z1Aw2_groups, &z0Aw1_list,
            &w0Az0_groups, &w1Az1_groups, &w1Az2_groups, &w0Az1_list,

            x, y, Ax, Ay, wgtx, wgty,
            x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
            y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list
        );
        
        /*--------------------------------------------------------------------*/
        // Alternate between neutral Upflip_v5s in z and strictly beneficial flips
        /*--------------------------------------------------------------------*/
        do
        {
            // Neutral Upflip_v5s
            Neutral_Upflip_v5s_In_z_And_Downflip_v5s_In_w1Az2_v5
            (
                z, Az, w, Aw, wgtz, wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                problem
            ) ;

            
            
            // Strictly beneficial flips
            FeasibleBeneficialFlips
            (
                x, y, Ax, Ay, wgtx, wgty, true, true,
                x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
                y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
                problem, &numFeasibleStrictlyBenFlips
            ) ;

        } while(numFeasibleStrictlyBenFlips) ;

        /*--------------------------------------------------------------------*/
        // Alternate between neutral dnflips in w and strictly beneficial flips
        // MOST LIKELY there are no neutral Downflip_v5s after strictly
        // beneficial flips, so no point for a while loop. Need to double check
        // though
        /*--------------------------------------------------------------------*/
        do
        {
            // Neutral Downflip_v5s
            Neutral_Downflips
            (
                w, Aw, z, Az, wgtw, wgtz,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                problem
            ) ;
            
            // Strictly beneficial flips
            FeasibleBeneficialFlips
            (
                x, y, Ax, Ay, wgtx, wgty, true, true,
                x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
                y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
                problem, &numFeasibleStrictlyBenFlips
            ) ;
        } while(numFeasibleStrictlyBenFlips) ;
    }
}

/*----------------------------------------------------------------------------*/

void Neutral_Upflip_v5s_In_z_And_Downflip_v5s_In_w1Az2_v5
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtz,
    VSPINT                       *wgtw,
    IdxGroupsByC         *z0Aw0_groups,
    Direct_Access_List     *z0Aw1_list,
    IdxGroupsByC         *z1Aw1_groups,
    IdxGroupsByC         *z1Aw2_groups,

    IdxGroupsByC         *w0Az0_groups,
    Direct_Access_List     *w0Az1_list,
    IdxGroupsByC         *w1Az1_groups,
    IdxGroupsByC         *w1Az2_groups,

    Problem               *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT indFlip ;
    char   varFlip ;

    VSPINT max_num_upflips = 0.01 * problem->graph->n ;

    VSPINT numFlips = 0 ;
    while( 1 )
    {
        if(numFlips >= max_num_upflips) break ;

        #ifdef DEBUG
            Check_All_Data
            (
                z, Az, *wgtz, w, Aw, *wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                problem
            ) ;
        #endif        

        // Choose a flip
        VSPINT randNum = (VSPINT) rand() ;
        Select_Random_Feasible_Neutral_Upflip_In_z_Or_w1Az2
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            randNum, problem, &indFlip, &varFlip
        ) ;

        #ifdef LOG_FLIPS_V5
        if( LOG_CONDITION )
        {
            if(indFlip != EMPTY)
            {
                char realVarFlip = (varFlip == 'z') ? z0Aw0_groups->var : w0Az0_groups->var ;
                char filename[150] ;
                snprintf(filename, sizeof(filename), "Flips_v5/problem_%s_Level_%lld.csv", problem->name, (long long) problem->level) ;
                Log(filename, "%c, %lld", realVarFlip, (long long) indFlip) ;
            }
        }
        #endif

        // Do the selected flip
        if(varFlip == 'z')
        {
            // Upflip_v5 in z
            VSP_ASSERT(z[indFlip] == 0) ;
            Upflip
            (
                z, Az, w, Aw, wgtz, wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

                indFlip, problem
            ) ;
            numFlips++ ;
        }
        else if (varFlip == 'w')
        {
            // Downflip_v5 in w
            VSP_ASSERT(w[indFlip] == 1) ;
            Downflip
            (
                w, Aw, z, Az, wgtw, wgtz,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            
                indFlip, problem
            ) ;
            numFlips++ ;
        }
        else
        {
            break ;
        }
        #ifdef DEBUG
            Check_All_Data
            (
                z, Az, *wgtz, w, Aw, *wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                problem
            ) ;
        #endif       
    }

    return ;
}



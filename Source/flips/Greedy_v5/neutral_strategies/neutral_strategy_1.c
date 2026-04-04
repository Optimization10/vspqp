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

#include "neutral_strategies_v5.h"


void Upflip
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

    VSPINT                     indFlip,
    Problem               *problem
) ;


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
        Assign_By_Weight_v5
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
            StrictlyBeneficialFlips_v5
            (
                x, y, Ax, Ay, wgtx, wgty,
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
            Neutral_Downflip_v5s_v5
            (
                w, Aw, z, Az, wgtw, wgtz,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                problem
            ) ;
            
            // Strictly beneficial flips
            StrictlyBeneficialFlips_v5
            (
                x, y, Ax, Ay, wgtx, wgty,
                x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
                y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
                problem, &numFeasibleStrictlyBenFlips
            ) ;
        } while(numFeasibleStrictlyBenFlips) ;
    }
}




/*----------------------------------------------------------------------------*/
// Do all feasible, neutral Downflip_v5s in z (first variable)
void Neutral_Downflip_v5s_v5
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
    SET_FUNCTION_COUNTER ;
    VSPINT gamma  = problem->graph->max_vtx_cost ;
    VSPINT *group =  z1Aw1_groups->groups[gamma] ;
    VSPINT *size  = &( z1Aw1_groups->sizes[gamma] ) ;
    VSPINT *c     = problem->graph->c ;
    VSPINT n      = problem->graph->n ;

    VSPINT numFlips = 0 ;
    while( 1 )
    { 
        if( *size )
        {
            VSPINT randNum = (VSPINT) rand() ;
            
            // Choose Downflip_v5 at random
            VSPINT indFlip = group[ randNum % *size ] ;

            #ifdef DEBUG
                // Check that the flip is neutral
                if( !IN_RANGE(indFlip, 0, n) ) 
                {
                    ERR("Flip out of bounds") ;
                }

                // Check that this is a neutral Downflip_v5
                VSP_ASSERT( z[indFlip]  == 1 ) ;
                VSP_ASSERT( Aw[indFlip] == 1 ) ;
                VSP_ASSERT( c[indFlip]  == gamma ) ;
            #endif

            #ifdef LOG_FLIPS_V5
                if( LOG_CONDITION )
                {
                    if(indFlip != EMPTY)
                    {
                        char filename[150] ;
                        snprintf(filename, sizeof(filename), "Flips_v5/problem_%s_Level_%lld.csv", problem->name, (long long) problem->level) ;
                        Log(filename, "%c, %lld", z1Aw1_groups->var, (long long) indFlip) ;
                    }
                }
            #endif

            // Do this flip
            Downflip_v5
            (
                z, Az, w, Aw, wgtz, wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

                indFlip, problem
            ) ;
            numFlips++ ;
        }
        else
        {
            #ifdef DEBUG
                // Check that there is no feasible, neutral Downflip_v5 in z
                for(VSPINT i = 0 ; i < n ; i++)
                {
                    if(z[i] == 1 && Aw[i] == 1 && c[i] == gamma)
                    {
                        ERR("Neutral Downflip_v5 exists but it was not detected") ;
                        break ;
                    }
                }
            #endif
            break ;
        }
    }
    return ;
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
            Check_All_Data_v5
            (
                z, Az, *wgtz, w, Aw, *wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                problem
            ) ;
        #endif        

        // Choose a flip
        VSPINT randNum = (VSPINT) rand() ;
        Select_Random_Feasible_Neutral_Upflip_v5_In_z_Or_w1Az2_v5
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
            Downflip_v5
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
            Check_All_Data_v5
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

void Select_Random_Feasible_Neutral_Upflip_v5_In_z_Or_w1Az2_v5
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

    VSPINT                     randNum,
    Problem               *problem,
    VSPINT                 *retIndFlip,
    char                   *retVarFlip
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT u = problem->u ;
    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    *retIndFlip = EMPTY ;
    *retVarFlip = '\0' ;


    // If wi = 1, Azi = 2 exists
    if(w1Az2_groups->minHeap->sizeHeap )
    {
        // Choose a random Downflip_v5 corresponding to the root. Most likely,
        // we will end up doing all of them, so randomness does not change
        // anything. But it's here for good measure
        VSPINT cval      = w1Az2_groups->minHeap->heap[1] ;
        VSPINT *group    = w1Az2_groups->groups[cval] ;
        VSPINT sizeGroup = w1Az2_groups->sizes[cval] ;

        // Don't do infeasible downflip
        *retIndFlip      = group[ randNum % sizeGroup ] ;
        *retVarFlip      = 'w' ;

        goto RET ;
    }

    // Select a feasible neutral Upflip_v5
    VSPINT sizeList = z0Aw1_list->size ;
    if(sizeList && *wgtz + gamma <= u)
    {
        *retIndFlip = z0Aw1_list->list[randNum % sizeList] ;
        *retVarFlip = 'z' ;
    }

RET :
    #ifdef DEBUG
        Select_Random_Feasible_Neutral_Upflip_v5_In_z_Or_w1Az2_v5_Check
        (
            z, Az, w, Aw, wgtz, wgtw, 
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups, 
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups, 
            randNum, problem, *retIndFlip, *retVarFlip
        );    
    #endif
    return ;
}




/*----------------------------------------------------------------------------*/
// Check:
// If there is a Downflip_v5 with wi = 1 and Azi >= 2, then the given flip is such
// a Downflip_v5
// If there is no such Downflip_v5, and there is a feasible, neutral Upflip_v5 in z,
// then the given flip is such a flip
// Else, if there is no such Upflip_v5 or Downflip_v5, the given EMPTY
/*----------------------------------------------------------------------------*/
void Select_Random_Feasible_Neutral_Upflip_v5_In_z_Or_w1Az2_v5_Check
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

    VSPINT                     randNum,
    Problem               *problem,
    VSPINT                givenIndFlip,
    char                  givenVarFlip
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;
    VSPINT u     = problem->u ;

    /*------------------------------------------------------------------------*/
    // Basic input validation
    /*------------------------------------------------------------------------*/
    if(!IN_RANGE(givenIndFlip, 0, n) && givenIndFlip != EMPTY)
    {
        ERR("Given flip out of bounds") ;
    }

    // Check feasibility
    if(givenVarFlip == 'z')
    {
        if(*wgtz + gamma > u)
        {
            ERR("Infeasible neutral Upflip_v5 selected") ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Actual checks
    /*------------------------------------------------------------------------*/

    bool is_w1Az2 = false ;
    VSPINT dnflipWitness = EMPTY ;

    bool is_z0Aw1cgamma = false ;
    VSPINT Upflip_v5Witness = EMPTY ;

    // If there is a w1Az2 flip but it's infeasible, we quit neutral
    // upflips early. Thus in this case, it's not an issue if there are still
    // neutral upflips left
    bool is_lb_infeasible_is_w1Az2 = false ;
    VSPINT lb_infeasible_witness = EMPTY ;

    // Search for any flips
    for(i = 0 ; i < n ; i++)
    {
        if(w[i] == 1 && Az[i] == 2 && *wgtw > c[i])
        {
            is_w1Az2      = true ;
            dnflipWitness = i ;
        }
        else if(w[i] == 1 && Az[i] == 2 && *wgtw == c[i])
        {
            is_lb_infeasible_is_w1Az2 = true ;
            lb_infeasible_witness = i ;
        }

        if(*wgtz + gamma <= u)
        {
            if(z[i] == 0 && Aw[i] == 1 && c[i] == gamma)
            {
                is_z0Aw1cgamma = true ;
                Upflip_v5Witness  = i ;
            }
        }
    }

    if(is_w1Az2)
    {
        if(givenVarFlip != 'w')
        {
            ERR("Downflip_v5 with wi = 1 and Azi = 2 not detected") ;
        }

        if(givenIndFlip != EMPTY)
        {
            if(w[givenIndFlip] != 1 || Az[givenIndFlip] != 2)
            {
                ERR("Malformed flip") ;
            }
        }
        else ERR("Given empty flip when there is one") ;
    }
    else if(is_z0Aw1cgamma)
    {
        // If there's a neutral upflips and there is no worry about making Az 3
        if(givenVarFlip != 'z' && is_lb_infeasible_is_w1Az2 == false)
        {
            ERR("Neutral flip not detected") ;
        }

        if(givenIndFlip != EMPTY)
        {
            if(z[givenIndFlip] != 0 || Aw[givenIndFlip] != 1 || c[givenIndFlip] != gamma)
            {
                ERR("Malformed flip") ;
            }
        }
        // else ERR("Given empty flip when there is one") ;
    }
    else
    {
        if(givenIndFlip != EMPTY || givenVarFlip != '\0')
        {
            ERR("No flip exists, but one is claimed") ;
        }
    }
    return ;
}
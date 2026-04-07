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
        Do neutral upflips in z. If S ever becomes nonempty, do the downflip(s)
        contained in S, and return to neutral upflips z. Keep doing this until
        there are no more neutral upflips in z and S is empty (see Rmk1 below)

        Do all strictly beneficial flips in (z, w), and let M be the number
        of strictly beneficial flips done
    End Repeat

    Initialize M = 1
    Repeat until M = 0:
        Do all neutral downflips in w

        Do all strictly beneficial flips in (z, w), and let M be the number of
        strictly beneficial flips done
    End Repeat
End Repeat
Do all strictly beneficial, and neutral downflips. Stop when there are no
strictly beneficial downflips available,

Rmk 1: It is theoretically possible that doing the downflip would violate the
       lower balance constraint. This is very unlikely, but I need to deal with
       it in the code. The code does not deal with this case at the moment. It
       will just do the flip and violate the lower balance constraint, which
       may cause other bugs down the line

Rmk 2: When we say "do all feasible, strictly beneficial flips", it is only
       the feasible, strictly beneficial flips visible to the data structures.
       It is possible that due to being close to the balance constraint, there
       ARE feasible, strictly beneficial upflips, but they do not correspond to
       the root of the max heap in the IdxGroupsByC data structs, so they are
       not done 
********************************************************************************

------------------------------------------------------------------------------*/
#include "core.h"
#include "greedy_flips_internal.h"
#include "utils.h"


/*----------------------------------------------------------------------------*/
// Check:
// If there is a downflip with wi = 1 and Azi >= 2, then the given flip is such
// a downflip
// If there is no such downflip, and there is a feasible, neutral upflip in z,
// then the given flip is such a flip
// Else, if there is no such upflip or downflip, the given EMPTY
/*----------------------------------------------------------------------------*/
void Select_Random_Feasible_Neutral_Upflip_In_z_Or_w1Az2_Check
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
            ERR("Infeasible neutral upflip selected") ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Actual checks
    /*------------------------------------------------------------------------*/

    bool is_w1Az2 = false ;
    UNUSED VSPINT dnflipWitness = EMPTY ;


    bool is_z0Aw1cgamma = false ;
    UNUSED VSPINT upflipWitness = EMPTY ;


    // If there is a w1Az2 flip but it's infeasible, we quit neutral
    // upflips early. Thus in this case, it's not an issue if there are still
    // neutral upflips left
    bool is_lb_infeasible_is_w1Az2 = false ;
    VSPINT lb_infeasible_witness = EMPTY ;

    // Search for any flips.
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
                upflipWitness  = i ;
            }
        }
    }

    if(is_w1Az2)
    {
        if(givenVarFlip != 'w')
        {
            ERR("Downflip with wi = 1 and Azi = 2 not detected") ;
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


void Select_Random_Feasible_Neutral_Upflip_In_z_Or_w1Az2
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
        Check_All_Data
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
        // Choose a random downflip corresponding to the root. Most likely,
        // we will end up doing all of them, so randomness does not change
        // anything. But it's here for good measure
        VSPINT cval      = w1Az2_groups->minHeap->heap[1] ;
        VSPINT *group    = w1Az2_groups->groups[cval] ;
        VSPINT sizeGroup = w1Az2_groups->sizes[cval] ;

        *retIndFlip      = group[ randNum % sizeGroup ] ;
        *retVarFlip      = 'w' ;

        goto RET ;
    }

    // Select a feasible neutral upflip
    VSPINT sizeList = z0Aw1_list->size ;
    if(sizeList && *wgtz + gamma <= u)
    {
        *retIndFlip = z0Aw1_list->list[randNum % sizeList] ;
        *retVarFlip = 'z' ;
    }

RET :
    #ifdef DEBUG
        Select_Random_Feasible_Neutral_Upflip_In_z_Or_w1Az2_Check
        (
            z, Az, w, Aw, wgtz, wgtw, 
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups, 
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups, 
            randNum, problem, *retIndFlip, *retVarFlip
        );    
    #endif
    return ;
}

void Neutral_Flips_Strategy_1
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT *wgtx,
    VSPINT *wgty,
    bool    deep_search,

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
        // Alternate between neutral upflips in z and strictly beneficial flips
        /*--------------------------------------------------------------------*/
        do
        {
            // Neutral upflips
            Neutral_Upflips_In_z_And_Downflips_In_w1Az2
            (
                z, Az, w, Aw, wgtz, wgtw,
                z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
                w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
                problem
            ) ;
            
            // Strictly beneficial flips, no deep search
            FeasibleBeneficialFlips
            (
                x, y, Ax, Ay, wgtx, wgty, deep_search, true,
                x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
                y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
                problem, &numFeasibleStrictlyBenFlips
            ) ;

        } while(numFeasibleStrictlyBenFlips) ;


        /*--------------------------------------------------------------------*/
        // Alternate between neutral dnflips in w and strictly beneficial flips
        // MOST LIKELY there are no neutral downflips after strictly
        // beneficial flips, so no point for a while loop. Need to double check
        // though
        /*--------------------------------------------------------------------*/
        do
        {
            // Neutral downflips
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
                x, y, Ax, Ay, wgtx, wgty, deep_search, true,
                x0Ay0_groups, x1Ay1_groups, x1Ay2_groups, x0Ay1_list,
                y0Ax0_groups, y1Ax1_groups, y1Ax2_groups, y0Ax1_list,
                problem, &numFeasibleStrictlyBenFlips
            ) ;
        } while(numFeasibleStrictlyBenFlips) ;
    }
}




/*----------------------------------------------------------------------------*/
// Do all feasible, neutral downflips in z (first variable)
// WARNING: I am not only doing neutral downflips here, but also strictly
// beneficial downflips
void Neutral_Downflips
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
            
            // Choose downflip at random
            VSPINT indFlip = group[ randNum % *size ] ;

            #ifdef DEBUG
                // Check that the flip is neutral
                if( !IN_RANGE(indFlip, 0, n) ) 
                {
                    ERR("Flip out of bounds") ;
                }

                // Check that this is a neutral downflip
                VSP_ASSERT( z[indFlip]  == 1 ) ;
                VSP_ASSERT( Aw[indFlip] == 1 ) ;
                VSP_ASSERT( c[indFlip]  == gamma ) ;
            #endif


            // Do this flip
            Downflip
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
                // Check that there is no feasible, neutral downflip in z
                for(VSPINT i = 0 ; i < n ; i++)
                {
                    if(z[i] == 1 && Aw[i] == 1 && c[i] == gamma)
                    {
                        ERR("Neutral downflip exists but it was not detected") ;
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

void Neutral_Upflips_In_z_And_Downflips_In_w1Az2
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

#ifdef DEBUG
        Check_All_Data
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
#endif 

    VSPINT indFlip ;
    char   varFlip ;

    VSPINT max_num_upflips = 0.01 * problem->graph->n ;

    VSPINT numFlips = 0 ;
    while( 1 )
    {


        if(numFlips >= max_num_upflips) break ;       

        // Choose a flip
        VSPINT randNum = (VSPINT) rand() ;
        Select_Random_Feasible_Neutral_Upflip_In_z_Or_w1Az2
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            randNum, problem, &indFlip, &varFlip
        ) ;



        // Do the selected flip
        if(varFlip == 'z')
        {
            // Upflip in z
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
            // Downflip in w
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
    }

    return ;
}










// Assign the lighter variable to z and the heavier one to w
void Assign_By_Weight
(
    VSPINT **z,
    VSPINT **w,
    VSPINT **Az,
    VSPINT **Aw,
    VSPINT **wgtz,
    VSPINT **wgtw,

    IdxGroupsByC **z0Aw0_groups,
    IdxGroupsByC **z1Aw1_groups,
    IdxGroupsByC **z1Aw2_groups,
    Direct_Access_List **z0Aw1_list,

    IdxGroupsByC **w0Az0_groups,
    IdxGroupsByC **w1Az1_groups,
    IdxGroupsByC **w1Az2_groups,
    Direct_Access_List **w0Az1_list,

    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT *wgtx,
    VSPINT *wgty,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list,

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list
)
{
    if (*wgtx < *wgty) {
        *z = x;
        *w = y;
        *Az = Ax;
        *Aw = Ay;
        *wgtz = wgtx;
        *wgtw = wgty;

        *z0Aw0_groups = x0Ay0_groups;
        *z1Aw1_groups = x1Ay1_groups;
        *z1Aw2_groups = x1Ay2_groups;
        *z0Aw1_list = x0Ay1_list;

        *w0Az0_groups = y0Ax0_groups;
        *w1Az1_groups = y1Ax1_groups;
        *w1Az2_groups = y1Ax2_groups;
        *w0Az1_list = y0Ax1_list;
    }
    else
    {
        *z = y;
        *w = x;
        *Az = Ay;
        *Aw = Ax;
        *wgtz = wgty;
        *wgtw = wgtx;

        *z0Aw0_groups = y0Ax0_groups;
        *z1Aw1_groups = y1Ax1_groups;
        *z1Aw2_groups = y1Ax2_groups;
        *z0Aw1_list = y0Ax1_list;

        *w0Az0_groups = x0Ay0_groups;
        *w1Az1_groups = x1Ay1_groups;
        *w1Az2_groups = x1Ay2_groups;
        *w0Az1_list = x0Ay1_list;
    }
}
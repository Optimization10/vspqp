#include "core.h"
#include "flips_private_v5.h"

#include "tmp.h"

void z0Aw1_List_Check_v5
(
    Direct_Access_List *List,
    VSPINT                *c,
    VSPINT                *z,
    VSPINT               *Aw,
    VSPINT             gamma,
    VSPINT                 n
)
{
    SET_FUNCTION_COUNTER ;
    WARN_IF_NOT_DEBUG_MODE ;

    #ifdef DEBUG
    Direct_Access_List_Check(List) ;
    VSPINT i, p ;
    VSPINT *list = List->list ;
    VSPINT size  = List->size ;
    
    // Get a witness array of all i in the list, and make sure all i satisfy
    // zi = 0, Awi = 1, ci = gamma
    VSPINT *witnessList = calloc(n, sizeof(VSPINT)) ;
    for(p = 0 ; p < size ; p++)
    {
        i = list[p] ;
        witnessList[i]++ ;
        if(z[i] != 0)
        {
            ERR("i with zi = 1 in z0Aw1 list") ;
        }

        if(Aw[i] != 1)
        {
            ERR("i with Awi != 1 in z0Aw1 list") ;
        }

        if(c[i] != gamma)
        {
            ERR("i with  ci != gamma in z0Aw1 list") ;
        }
    }

    VSPINT *witnessAccurate = calloc(n, sizeof(VSPINT)) ;
    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == 0 && Aw[i] == 1 && c[i] == gamma)
        {
            witnessAccurate[i]++ ;
        }
    }

    for(i = 0 ; i < n ; i++)
    {
        if(witnessAccurate[i] != witnessList[i])
        {
            ERR("Mismatch between witnesses") ;
        }

        if(witnessList[i] != 0 && witnessList[i] != 1)
        {
            ERR("An i is in a list more than once") ;
        }

        if(witnessAccurate[i] != 0 && witnessAccurate[i] != 1)
        {
            ERR("Very strange error") ;
        }
    }


    free(witnessList) ;
    free(witnessAccurate) ;
    #endif
    return ;
}




// In addition to checking the data,check that there are no strictly beneficial
// flips 
void Check_Outputs_v5
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    Problem       *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;
    Check_Inputs(x, Ax, wgtx, y, Ay, wgty, problem) ;

    VSPINT n = problem->graph->n ;
    VSPINT u = problem->u ;
    VSPINT *c = problem->graph->c ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    VSPINT ps[4], cs[4], witness[4], isFeas[4], gains[4] ;
    Compute_ps_And_cs_Dbg_v5(x, Ax, wgtx, y, Ay, wgty, problem, ps, cs, witness) ;

    Compute_Feasibility_And_Gains_v5
    (
        Ax, Ay, wgtx, wgty, problem,
        witness, isFeas, gains
    ) ;

    for(VSPINT k = 0 ; k < 4 ; k++)
    {
        if(isFeas[k] && gains[k] > 0)
        {
            ERR("Quit strictly beneficial flips with one still available") ;
        }
    }

    // WARNING: this may give false positive, as we exit when structs tell
    // us there are no more flips, but remember non-root heap node may
    // correspond to a beneficial Upflip_v5 that the structs can't see
    // for(VSPINT i = 0 ; i < n ; i++)
    // {
    //     if(x[i])
    //     {
    //         if(GAIN_FLIP(-1, c[i], gamma, Ay[i]) > 0 && FEASIBLE_DNFLIP(wgtx, c[i]))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Downflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    //     else
    //     {
    //         if(GAIN_FLIP(1, c[i], gamma, Ay[i]) > 0 && FEASIBLE_UPFLIP(wgtx, c[i], u))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Upflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }

    //     if(y[i])
    //     {
    //         if(GAIN_FLIP(-1, c[i], gamma, Ax[i]) > 0 && FEASIBLE_DNFLIP(wgty, c[i]))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Downflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    //     else
    //     {
    //         if(GAIN_FLIP(1, c[i], gamma, Ax[i]) > 0 && FEASIBLE_UPFLIP(wgty, c[i], u))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Upflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    // }
}





#include "core.h"
/*----------------------------------------------------------------------------*/
// This file is for functions which are commonly used throughout many of
// the flip codes, whether for debugging or no. I will try to refrain from
// data structures as much as possible here in order to maintain flexibility
/*----------------------------------------------------------------------------*/

#include "shared_api.h"
#include "flips_common_api.h"




void Compute_ps_And_cs_Dbg
(
    VSPINT                    *x,
    VSPINT                   *Ax,
    VSPINT                  wgtx,
    VSPINT              *xlocked, // Ignored if NULL
    VSPINT                    *y,
    VSPINT                   *Ay,
    VSPINT                  wgty,
    VSPINT              *ylocked, // Ignored if NULL
    VSP_problem         *problem,
    VSPINT                 ps[4],
    VSPINT                 cs[4],
    VSPINT          witnesses[4]
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    // Check that either one or both variables are locked
    if((xlocked == NULL && ylocked != NULL) ||
       (xlocked != NULL && ylocked == NULL)    )
    {
        ERR("Only one variable locked") ;
    }

    VSPINT i, xi, yi, Axi, Ayi, ci ;
    VSPINT n     = problem->Graph->n ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;
    VSPINT *c    = problem->Graph->c ;
    VSPINT  u    = problem->u ;

    /*------------------------------------------------------------------------*/
    // Get the p values
    /*------------------------------------------------------------------------*/
    VSPINT px0 = n + 1, py0 = n + 1, px1 = -1, py1 = -1 ;

    for(i = 0 ; i < n ; i++)
    {
        xi  = x[i] ;
        yi  = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;

        if( UNLOCKED(xlocked, i) )
        {
            // Record px's
            if(xi)
            {
                if(Ayi > px1)
                {
                    px1 = Ayi ;
                }
            }
            else
            {
                if(Ayi < px0)
                {
                    px0 = Ayi ;
                }
            }
        }

        if( UNLOCKED(ylocked, i) )
        {
            // Record py's
            if(yi)
            {
                if(Axi > py1)
                {
                    py1 = Axi ;
                }
            }
            else
            {
                if(Axi < py0)
                {
                    py0 = Axi ;
                }
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Record the c's and an index for each c
    // Also record:
    // x0witness, y0witness, x1witness, y1witness - indices achieving the
    // optimal p and c values
    /*------------------------------------------------------------------------*/
    VSPINT cx0 = 0, cy0 = 0, cx1 = gamma + 1, cy1 = gamma + 1 ;
    VSPINT x0witness = EMPTY, y0witness = EMPTY ;
    VSPINT x1witness = EMPTY, y1witness = EMPTY ;

    for(i = 0 ; i < n ; i++)
    {
        xi = x[i] ;
        yi = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;
        ci = c[i] ;

        // Deal with x
        if(xi == 0 && Ayi == px0 && UNLOCKED(xlocked, i) )
        {
            if(ci > cx0 )
            {
                cx0 = ci ;
                x0witness = i ;
            }
        }

        if(xi == 1 && Ayi == px1 && UNLOCKED(xlocked, i))
        {
            if(ci < cx1 )
            {
                cx1 = ci ;
                x1witness = i ;
            }
        }

        // Deal with y
        if(yi == 0 && Axi == py0 && UNLOCKED(ylocked, i))
        {
            if(ci > cy0 )
            {
                cy0 = ci ;
                y0witness = i ;
            }
        }

        if(yi == 1 && Axi == py1 && UNLOCKED(ylocked, i) )
        {
            if(ci < cy1 )
            {
                cy1 = ci ;
                y1witness = i ;
            }
        }
    }

    // Prepare return
    ps[X0] = px0 ;
    ps[X1] = px1 ;
    ps[Y0] = py0 ;
    ps[Y1] = py1 ;

    cs[X0] = cx0 ;
    cs[X1] = cx1 ;
    cs[Y0] = cy0 ;
    cs[Y1] = cy1 ;

    witnesses[X0] = x0witness ;
    witnesses[X1] = x1witness ;
    witnesses[Y0] = y0witness ;
    witnesses[Y1] = y1witness ;

    return ;
}




/*------------------------------------------------------------------------------
Record the gain even if it's negative

inds[k] - k is for the encoding, inds[k] is a vertex number. inds[k] may be
          EMPTY

retIsFeas[k] - if inds[k] != EMPTY and the flip is feasible, the it's 1
               otherwise, it's 0

retGains[k] - if retIsFeas[k] = 1, then this is the gain of the flip. If
              retIsFeas[k] = 0, then this is garbage
------------------------------------------------------------------------------*/
void Compute_Feasibility_And_Gains
(
    VSPINT                *Ax,
    VSPINT                *Ay,
    VSPINT               wgtx,
    VSPINT               wgty,
    VSP_problem      *problem,
    VSPINT             inds[4],
    VSPINT        retIsFeas[4],
    VSPINT         retGains[4]
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i ;
    VSPINT n  = problem->Graph->n ;
    VSPINT *c = problem->Graph->c ;
    VSPINT u  = problem->u ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;

    for(i = 0 ; i < 4 ; i++) retIsFeas[i] = 0 ;

    if(inds[X0] != EMPTY && FEASIBLE_UPFLIP(wgtx, c[inds[X0]], u) )
    {
        retIsFeas[X0] = 1 ;
        retGains [X0] = GAIN_FLIP(1,  c[inds[X0]], gamma, Ay[inds[X0]]) ;
    }

    if(inds[Y0] != EMPTY && FEASIBLE_UPFLIP(wgty, c[inds[Y0]], u) )
    {
        retIsFeas[Y0] = 1 ;
        retGains [Y0] = GAIN_FLIP(1,  c[inds[Y0]], gamma, Ax[inds[Y0]]) ;
    }

    if(inds[X1] != EMPTY && FEASIBLE_DNFLIP(wgtx, c[inds[X1]]) )
    {
        retIsFeas[X1] = 1 ;
        retGains [X1] = GAIN_FLIP(-1,  c[inds[X1]], gamma, Ay[inds[X1]]) ;
    }

    if(inds[Y1] != EMPTY && FEASIBLE_DNFLIP(wgty, c[inds[Y1]]) )
    {
        retIsFeas[Y1] = 1 ;
        retGains [Y1] = GAIN_FLIP(-1,  c[inds[Y1]], gamma, Ax[inds[Y1]]) ;
    }

    return ;
}


void Do_Flip_Basic
(
    VSPINT    *z,
    VSPINT   *Az,
    VSPINT     i,
    VSPINT *wgtz,
    VSPINT    *c,
    VSPINT   *Ai,
    VSPLONG *Ap
)
{
    VSPINT dir, j ;
    VSPLONG end, p ;

    // Do the flip and update weight
    dir    = 1 - 2 * z[i] ;
    z[i]  += dir ;
    *wgtz += dir * c[i] ;

    // Update Az
    end   = Ap[i+1] ;
    for(p = Ap[i] ; p < end ; p++)
    {
        j = Ai[p] ;
        Az[ j ] += dir ;
    }
    Az[i] += dir ;
}
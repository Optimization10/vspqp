#ifndef FLIPS_COMMON_API_H
#include "shared_api.h"

typedef struct
{
    char   var ;
    VSPINT i;
    VSPINT gain ;
    VSPINT ci ;
    VSPINT zi ;
    VSPINT Awi ;
    VSPINT gamma ;
    VSPINT wgt ;
    VSPINT u ;
    VSPINT dir ;
} FlipDbg ;

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
) ;

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
) ;

void GeneralFlipSelector_Check
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                    *xlocked, // Ignored if NULL
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    VSPINT                    *ylocked, // Ignored if NULL
    VSPINT               givenFlipType, // NONE, X0, Y0, X1, Y1
    VSPINT                givenIndFlip, 
    VSP_problem               *problem,
    Obj_Effect           expFlipEffect
) ;

void Do_Flip_Basic
(
    VSPINT    *z,
    VSPINT   *Az,
    VSPINT     i,
    VSPINT *wgtz,
    VSPINT    *c,
    VSPINT   *Ai,
    VSPINT   *Ap
) ;

#endif
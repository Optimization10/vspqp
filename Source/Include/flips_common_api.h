#ifndef FLIPS_COMMON_API_H
#include "shared_api.h"

#define NO_FLIP -1
#define XUP 0
#define YUP 1
#define XDN 2
#define YDN 3

#define NONE -1
#define X0 0
#define Y0 1
#define X1 2
#define Y1 3

#define NO_VAR -1
#define X_VAR   0
#define Y_VAR   1

// Plug in 0, 1, 2, or 3 according to our encoding
// X0 --->  1
// Y0 --->  1
// X1 ---> -1
// Y1 ---> -1
#define GET_DIR(x) ( ( (x) < 2 ) + (x >= 2)*(-1) )

// Plug in 0, 1, 2, or 3 according to our encoding
// X0 ---> 0
// Y0 ---> 0
// X1 ---> 1
// Y1 ---> 1
#define GET_ZI(x) ( x >= 2 )

// X0 = 0 --> 1 = Y_VAR
// Y0 = 1 --> 0 = X_VAR
// X1 = 2 --> 1 = Y_VAR
// Y1 = 3 --> 0 = X_VAR
#define GET_OTHER_VAR(flip)   ( ( (flip) + 1) % 2 )
#define IS_UPFLIP(flip) ( (flip) <= 1 )
#define IS_DNFLIP(flip) ( (flip) >= 2 )

// 0 --> 1
// 1 --> 0
#define OTHER_VAR(v) ( 1 - (v) )

// Here, x is one of X0, Y0, X1, Y1
// X0, X1 (that is 0 and 2) map to 0 (X_VAR), and Y0 and Y1 (that is 1 and 3)
// map to 1 (Y_VAR) 
#define GET_VAR(x) ( (x) % 2)

#define UNLOCKED(locked, i) (  ((locked) == NULL) || ( !(locked)[(i)] )  )

#define FEASIBLE_UPFLIP(wgt, ci, u) ( (wgt) + (ci) <= (u) )

// Lower bound of 1 is disabled
#define FEASIBLE_DNFLIP(wgt, ci) ( 1 )

#define IS_FEASIBLE_FLIP(zi, wgt, ci, u) ( zi ? FEASIBLE_DNFLIP(wgt, ci) : FEASIBLE_UPFLIP(wgt, ci, u) )


// zi = 0 -->  1
// zi = 1 --> -1
#define DIR_FLIP(zi) ( 1 - 2 * (zi) )

#define GAIN_FLIP(dir, ci, gamma, Awi) ( (dir) * ( (ci) - (gamma) * (Awi) ) )

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
    Problem             *problem,
    VSPINT                 ps[4],
    VSPINT                 cs[4],
    VSPINT          witnesses[4]
) ;

void Compute_Feasibility_And_Gains
(
    VSPINT                 *Ax,
    VSPINT                 *Ay,
    VSPINT                wgtx,
    VSPINT                wgty,
    Problem           *problem,
    VSPINT             inds[4],
    VSPINT        retIsFeas[4],
    VSPINT         retGains[4]
) ;

void Do_Flip_Basic
(
    VSPINT    *z,
    VSPINT   *Az,
    VSPINT     i,
    VSPINT *wgtz,
    VSPINT    *c,
    VSPINT   *Ai,
    VSPLONG  *Ap
) ;

#endif
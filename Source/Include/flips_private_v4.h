#ifndef FLIPS_PRIVATE_V4_H
#define FLIPS_PRIVATE_V4_H
#include "shared_api.h"
#include "ds_api.h"
#include "flips_public_v4.h"
#include "flips_common_api.h"
#define LOG_CONDITION (problem->level == 0)


/*----------------------------------------------------------------------------*/
/* Miscellaneous                                                              */
/*----------------------------------------------------------------------------*/


typedef enum
{
    IMPROVE_OBJ_ONLY,
    IMPROVE_OBJ_OR_NEUTRAL_DOWNFLIP,
    GAIN_GEQ_0
} Expected_Flip_Effect ;

void Initialize_All_Structs
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    IdxGroupsByC         *x0Ay0_groups,
    Direct_Access_List     *x0Ay1_list,
    IdxGroupsByC         *x1Ay1_groups,
    IdxGroupsByC         *x1Ay2_groups,

    IdxGroupsByC         *y0Ax0_groups,
    Direct_Access_List     *y0Ax1_list,
    IdxGroupsByC         *y1Ax1_groups,
    IdxGroupsByC         *y1Ax2_groups,

    Problem               *problem
) ;

void SelectRandomOptimalStrictlyBeneficialFlip
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list, // only used in debugging

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list, // only used in debugging

    Problem       *problem,
    char           *retFlipVar, // '\0' if no flip
    VSPINT         *retIndFlip  // EMPTY if no flip
) ;

void Do_Flip
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                       *wgtz,
    VSPINT                          *w,
    VSPINT                         *Aw,
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

void Downflip
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtz,
    VSPINT                       *wgtw, // used for debugging
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

// Assume that upflip happens in z
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

    Problem                 *problem,
    VSPINT              *retNumFlipsDone
) ;

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

    Problem                 *problem,
    VSPINT              *retNumFlipsDone
) ;

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
) ;

void Neutral_Flips_Strategy_1
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
) ;





/*----------------------------------------------------------------------------*/
/* Flip selectors                                                             */
/*----------------------------------------------------------------------------*/
void SelectRandomOptimalBeneficialFlip
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list, // only used in debugging

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list, // only used in debugging

    Problem       *problem,
    char           *retFlipVar, // '\0' if no flip
    VSPINT         *retIndFlip  // EMPTY if no flip
) ;


/*----------------------------------------------------------------------------*/
/* Debugging                                                                  */
/*----------------------------------------------------------------------------*/
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
} Flip ;




void z0Aw1_List_Check
(
    Direct_Access_List *List,
    VSPINT                *c,
    VSPINT                *z,
    VSPINT               *Aw,
    VSPINT             gamma,
    VSPINT                 n
) ;

void FeasibleStrictlyBeneficialFlipSelector_Check
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    char                  givenVarFlip,
    VSPINT                givenIndFlip,
    Problem               *problem,
    Expected_Flip_Effect expFlipEffect
) ;

void Check_All_Data
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    IdxGroupsByC         *x0Ay0_groups,
    Direct_Access_List     *x0Ay1_list,
    IdxGroupsByC         *x1Ay1_groups,
    IdxGroupsByC         *x1Ay2_groups,

    IdxGroupsByC         *y0Ax0_groups,
    Direct_Access_List     *y0Ax1_list,
    IdxGroupsByC         *y1Ax1_groups,
    IdxGroupsByC         *y1Ax2_groups,

    Problem               *problem
) ;


// ASSUMED: The flip has NOT happened yet
// The "best feasible flip according to structs" is a flip which achieves the
// gain (accounting for feasibility)
// max{ c_(x0), c_(y0), gamma * p_(x1) - c_(x1), gamma * p_(y1) - c_(y1) }
// This function checks that the given flup achieves this gain
// If there is no flip with feasible gain according to the structs, then
// the function expects givenVarFlip = '\0' and givenIndFlip = EMPTY 
void FeasibleBeneficialFlipSelector_Check
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    char          givenVarFlip,
    VSPINT        givenIndFlip,
    Problem       *problem
) ;

// Designed to be called at the very beginning of flips
void Check_Inputs
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    Problem       *problem
) ;

// Designed to be called at the very end of flips
// In addition to checking the data,check that there are no strictly beneficial
// flips 
void Check_Outputs
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    Problem       *problem
) ;

void Compute_Max_Group_Sizes
(
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT gamma,
    VSPINT n
) ;

// Let z be the lighter variable and w the heavier one
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
) ;

/*------------------------------------------------------------------------------
Encoding: 
0 --> x0
1 --> y0
2 --> x1
3 --> y1

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
    Problem      *problem,
    VSPINT             inds[4],
    VSPINT        retIsFeas[4],
    VSPINT         retGains[4]
) ;



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
) ;

#endif
#include "ds_api.h"

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
    Problem                   *problem
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
    Problem                   *problem
) ;


void Do_Flip
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                       *wgtz,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtw, // used for debugging only
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

void Compute_Max_Group_Sizes
(
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT gamma,
    VSPINT n
) ;

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
) ;

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
) ;


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


void SelectRandomOptimalBeneficialFlip
(
    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,
    bool deep_search,
    bool strict,

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


void FeasibleBeneficialFlips
(
    VSPINT                            *x,
    VSPINT                            *y,
    VSPINT                           *Ax,
    VSPINT                           *Ay,
    VSPINT                         *wgtx,
    VSPINT                         *wgty,
    bool                     deep_search,
    bool                          strict,

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
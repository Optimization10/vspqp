#ifndef FM_FLIPS_INTERNAL_H
#define FM_FLIPS_INTERNAL_H

#include "shared_api.h"
#include "FM_flips_api.h"
#include "ds_api.h"
#include "flips_api.h"
#include "flips_common_api.h"

#define ZDN 0
#define WDN 1
#define ZUP 2
#define WUP 3

#define PRINT_FM_FLIPS

#define GRADI(ci, gamma, Awi) ( ((ci) - (gamma) * (Awi)) )


typedef struct
{
    VSPINT *z     ;
    VSPINT *Az    ;
    VSPINT *Aw    ;
    VSPINT *c     ;
    VSPINT  n     ;
    VSPINT  gamma ;
    VSPINT weight ;
    VSPINT *Ai ;
    VSPINT *Ap ;
    VSPINT u ;

    VSPINT p_lo ; // min (Aw)_i s.t. z_i = 0 and i is unlocked
    VSPINT p_hi ; // max (Aw)_i s.t. z_i = 1 and i is unlocked

    VSPINT c_hi ; // max c_i s.t. z_i = 0, (Aw)_i = p_lo, and i is unlocked
    VSPINT c_lo ; // min c_i s.t. z_i = 1, (Aw)_i = p_hi, and i is unlocked

    Doubly_LLs_Unsized *L0 ;
    Doubly_LLs_Unsized *L1 ;

    Direct_Access_List *O_up ; // {i : z_i = 0, (Aw)_i = p_lo, c_i = c_hi, i unlocked}
    Direct_Access_List *O_dn ; // {i : z_i = 1, (Aw)_i = p_hi, c_i = c_lo, i unlocked}

    Max_Queue_v1 *A_up ; // {i : z_i = 0, (Aw)_i = p_lo, c_i < c_hi, i unlocked}
    Min_Queue_v1 *A_dn ; // {i : z_i = 1, (Aw)_i = p_hi, c_i > c_lo, i unlocked}

    VSPINT *locked ; // 0 if unlocked, 1 if locked

    char var ; // z or w

} FM_Flips_Struct ;


/*----------------------------------------------------------------------------*/
/*      Update When zi from 0 to 1                                            */
/*----------------------------------------------------------------------------*/
void FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_up
) ;

void FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_dn
) ;

FM_Flips_Struct *Construct_FM_Flips_Struct
(
    VSPINT            *z,
    VSPINT           *Az,
    VSPINT           *Aw,
    VSPINT            *c,
    VSPINT             n,
    VSPINT        weight,
    VSPINT         gamma,
    VSPINT           *Ai,
    VSPINT           *Ap,
    VSPINT             u
) ;

void Destruct_FM_Flips_Struct
(
    FM_Flips_Struct *data
) ;

void FM_Flips_Struct_Build
(
    FM_Flips_Struct *data,
    VSPINT         gamma
) ;

void Transfer_Nodes_From_Nonempty_A_up_To_The_Empty_O_Up
(
    Max_Queue_v1        *A_up,
    Direct_Access_List  *O_up
) ;


// Take the linked list corresponding to which_ll and builds O_up, A_up, and c_hi
// The function resets O_up and A_up so it is not neccesary to reset them before
void Build_O_up_And_A_up_And_c_hi_From_Linked_List
(
    FM_Flips_Struct    *data,
    VSPINT          which_ll
) ;

void Build_O_dn_And_A_dn_And_c_lo_From_Linked_List
(
    FM_Flips_Struct    *data,
    VSPINT          which_ll
) ;

void Set_New_p_lo_And_c_hi_And_O_up_And_A_up_When_zi_0_To_1
(
    FM_Flips_Struct        *data,
    VSPINT              old_p_lo
) ;

void FM_Flips_Struct_Check
(
    FM_Flips_Struct *data,
    bool         check_Az
) ;

void FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void Transfer_Nodes_From_Nonempty_A_dn_To_The_Empty_O_dn
(
    Min_Queue_v1        *A_dn,
    Direct_Access_List  *O_dn
) ;

void FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

Flip FM_Flips_Do_Best_Flip
(
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;

void FM_Flips_Check_Best_Flip
(
    VSPINT              *z,
    VSPINT              *w,
    VSPINT             *Az,
    VSPINT             *Aw,
    VSPINT       *locked_z,
    VSPINT       *locked_w,
    VSP_problem   *problem,
    VSPINT        ind_flip, // EMPTY if none
    char          var_flip
) ;

/*----------------------------------------------------------------------------*/
// Version 2
Flip FM_Flips_Do_Best_Flip_v2
(
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;


void FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v2
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_up
) ;

void FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v2
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_dn
) ;

void FM_Flips_Struct_Reset
(
    FM_Flips_Struct *data
) ;

void FM_Flips_Struct_Check_v2
(
    FM_Flips_Struct *data,
    bool         check_Az
) ;

void FM_Flips_Check_Best_Flip_v2
(
    VSPINT                *z,
    VSPINT                *w,
    VSPINT               *Az,
    VSPINT               *Aw,
    VSPINT         *locked_z,
    VSPINT         *locked_w,
    VSP_problem     *problem,
    Flip          flip_given
) ;

VSPINT SelectBestFlipRandomly
(
    VSPINT numFlips[4] 
) ;

/*----------------------------------------------------------------------------*/
// VERSION 3
Flip FM_Flips_Do_Best_Flip_v3
(
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;


void FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v3
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v3
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v3
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v3
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v3
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_up
) ;

void FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v3
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_dn
) ;

void FM_Flips_Struct_Reset
(
    FM_Flips_Struct *data
) ;

void FM_Flips_Struct_Check_v3
(
    FM_Flips_Struct *data,
    bool         check_Az
) ;

void FM_Flips_Check_Best_Flip_v3
(
    VSPINT                *z,
    VSPINT                *w,
    VSPINT               *Az,
    VSPINT               *Aw,
    VSPINT         *locked_z,
    VSPINT         *locked_w,
    VSP_problem     *problem,
    Flip          flip_given
) ;

void FM_Flips_Roll_Back_v3
(
    VSP_state                 *Z,
    VSP_state                 *W,
    FM_Flips_Struct    *Z_struct,
    FM_Flips_Struct    *W_struct,
    Flip                 *memory,
    VSP_problem         *problem,
    VSP_options         *options,
    VSPINT                nflips,
    VSPINT               locStop  // undo flips NOT including locStop
) ;

void FM_Flips_Do_Best_Flip_According_To_DS_Until_Fixed_Dry_Run_And_Roll_Back_To_Best_Location_v3
(
    VSP_state                 *Z,
    VSP_state                 *W,
    FM_Flips_Struct    *Z_struct,
    FM_Flips_Struct    *W_struct,
    VSPINT          lengthDryRun,
    VSP_problem         *problem,
    VSP_options         *options
) ;



VSPINT Flips_SelectBestFlipRandomly_v2
(
    VSPINT numFlips[4] 
) ;

int SelectBestFlipRandomlyBiased
(
    VSPINT numFlips[4],
    VSPINT        wgtz,
    VSPINT        wgtw 
) ;

void PrintFlipInfo
(
    VSPINT    runningGain,
    Flip             flip,
    int           dirFlip, // 1 for upflip, -1 for downflip
    VSP_problem  *problem
) ;

void WriteToFMFlipsLog
(
    char *msg
) ;



/*----------------------------------------------------------------------------*/
// Version 6
Flip FM_Flips_Do_Best_Flip_v6
(
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;


void FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v6
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v6
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v6
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;


void FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v6
(
    FM_Flips_Struct *data,
    VSPINT              j,
    VSPINT        Azj_old
) ;

void FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v6
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_up
) ;

void FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v6
(
    FM_Flips_Struct *data,
    VSPINT              i,
    VSPINT          idx_i  // index of i in O_dn
) ;

void FM_Flips_Struct_Reset
(
    FM_Flips_Struct *data
) ;

void FM_Flips_Struct_Check_v6
(
    FM_Flips_Struct *data,
    bool         check_Az
) ;

void FM_Flips_Check_Best_Flip_v6
(
    VSPINT                *z,
    VSPINT                *w,
    VSPINT               *Az,
    VSPINT               *Aw,
    VSPINT         *locked_z,
    VSPINT         *locked_w,
    VSP_problem     *problem,
    Flip          flip_given
) ;

VSPINT SelectBestFlipRandomly
(
    VSPINT numFlips[4] 
) ;

#endif
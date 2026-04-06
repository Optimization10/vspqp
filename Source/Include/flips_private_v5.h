#ifndef FLIPS_PRIVATE_V5_H
#define FLIPS_PRIVATE_V5_H
#include "shared_api.h"
#include "ds_api.h"
#include "flips_public_v5.h"
#define INSERT_IN_LIST(newListMem, insertLoc, list, pointers) \
                        (list)[ (insertLoc) ] = (newListMem) ; \
                        (pointers)[ (newListMem) ] = (insertLoc)

#define LOG_CONDITION (problem->level == 0)

/*----------------------------------------------------------------------------*/
/* Heaps                                                                      */
/*----------------------------------------------------------------------------*/
// Heap containing unique values of c
// Expects entries in 1, 2, ..., gamma
typedef struct
{
    VSPINT *heap ;
    VSPINT *heapPtrs ; // -1 if not in the heap
    VSPINT  sizeHeap ;
    VSPINT  gamma ; // Max size of the heap
    VSPINT *c ;
}  Min_Heap_v5 ;

// Heap containing unique values of c
// Expects entries in 1, 2, ..., gamma
typedef struct
{   
    VSPINT *heap ;
    VSPINT *heapPtrs ;
    VSPINT  sizeHeap ;
    VSPINT  gamma ; // Max size of the heap
    VSPINT *c ;
} Max_Heap_v5 ;

// Min Heap functions
Min_Heap_v5 *Min_Heap_Init_v5(VSPINT gamma, VSPINT *c);
void Min_Heap_Destroy_v5(Min_Heap_v5 *heap);
void Min_Heap_Check_v5(Min_Heap_v5 *Heap);
void Min_Heap_Bubble_Down_v5(Min_Heap_v5 *Heap, VSPINT p);
void Min_Heap_Bubble_Up_v5(Min_Heap_v5 *Heap, VSPINT p, VSPINT valToBubble);
void Min_Heap_Build_v5(Min_Heap_v5 *Heap);
void Min_Heap_Insert_v5(Min_Heap_v5 *Heap, VSPINT valToInsert);
void Min_Heap_Delete_v5(Min_Heap_v5 *Heap, VSPINT toDel);

// Max Heap functions
Max_Heap_v5 *Max_Heap_Init_v5(VSPINT gamma, VSPINT *c);
void Max_Heap_Destroy_v5(Max_Heap_v5 *heap);
void Max_Heap_Check_v5(Max_Heap_v5 *Heap);
void Max_Heap_Bubble_Down_v5(Max_Heap_v5 *Heap, VSPINT p);
void Max_Heap_Bubble_Up_v5(Max_Heap_v5 *Heap, VSPINT p, VSPINT valToBubble);
void Max_Heap_Build_v5(Max_Heap_v5 *Heap);
void Max_Heap_Insert_v5(Max_Heap_v5 *Heap, VSPINT valToInsert);
void Max_Heap_Delete_v5(Max_Heap_v5 *Heap, VSPINT toDel);


// /*----------------------------------------------------------------------------*/
// /* IdxGroupsByC                                                               */
// /*----------------------------------------------------------------------------*/
// typedef enum
// {
//     MIN_HEAP_FLAG,
//     MAX_HEAP_FLAG
// } MinMax_Heap_Flag ;

// typedef struct
// {
//     Min_Heap_v5 *minHeap;
//     Max_Heap_v5 *maxHeap;
//     VSPINT      **groups ;
//     VSPINT      *sizes ;
//     VSPINT      *maxSizes ; // kth entry is num indices in c with ci = k
//     VSPINT      *ptrs ;     // -1 for entries which are not here
//     VSPINT       gamma ;

//     const VSPINT *c ;
//     const VSPINT *z ;
//     const VSPINT *Aw ;
//     VSPINT targetzi ;
//     VSPINT targetAwi ;
//     char   var ;

//     MinMax_Heap_Flag flag ; // -1 for min, 1 for max
// } IdxGroupsByC ;

IdxGroupsByC *IdxGroupsByC_Init
(
    VSPINT  gamma,
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT nvtx, // Number of vertices at current level
    VSPINT targetzi,
    VSPINT targetAwi,
    VSPINT *z,
    VSPINT *Aw,
    char   var,
    MinMax_Heap_Flag flag
) ;

// void IdxGroupsByC_Destroy(IdxGroupsByC *Groups) ;

void IdxGroupsByC_Insert
(
    IdxGroupsByC *Groups,
    VSPINT idxToInsert
) ;

void IdxGroupsByC_Delete
(
    IdxGroupsByC *Groups,
    VSPINT idxToDel
) ;

void IdxGroupsByC_See_Top_Cval_If_Feasible
(
    IdxGroupsByC *Groups,
    VSPINT wgt, // Current weight
    VSPINT u,
    VSPINT *cvalRet  // EMPTY if not feasible
) ;

void IdxGroupsByC_Get_Random_Idx_Top_C_v5
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT       *idxReturn  // Return
) ;

void IdxGroupsByC_See_Random_Idx_Top_C
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT       *idxReturn  // Return
) ;

void IdxGroupsByC_Check
(
    IdxGroupsByC *Groups,
    VSPINT             n  // Size of z and Aw arrays
) ;

void IdxGroupsByC_Num_Idx_With_Top_C
(
    IdxGroupsByC    *Groups,
    VSPINT       *retNumIdx
) ;

void IdxGroupsByC_See_Top_Feasible_Cval_v5
(
    IdxGroupsByC      *z0Aw0_groups,
    Problem            *problem,
    VSPINT                     wgtz,
    VSPINT          *retMaxFeasCval,
    VSPINT                       *z, // Used for debugging only
    VSPINT                      *Aw // Used for debugging only
) ;

/*----------------------------------------------------------------------------*/
/* Miscellaneous                                                              */
/*----------------------------------------------------------------------------*/


typedef enum
{
    IMPROVE_OBJ_ONLY,
    IMPROVE_OBJ_OR_NEUTRAL_Downflip_v5,
    GAIN_GEQ_0
} Expected_Flip_Effect_v5 ;

void Initialize_All_Structs_v5
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

void SelectRandomOptimalStrictlyBeneficialFlip_v5
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

void Do_Flip_v5
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


void StrictlyBeneficialFlips_v5
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

void FeasibleBeneficialFlips_v5
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
) ;





/*----------------------------------------------------------------------------*/
/* Flip selectors                                                             */
/*----------------------------------------------------------------------------*/
void SelectRandomOptimalBeneficialFlip_v5
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


void z0Aw1_List_Check_v5
(
    Direct_Access_List *List,
    VSPINT                *c,
    VSPINT                *z,
    VSPINT               *Aw,
    VSPINT             gamma,
    VSPINT                 n
) ;

void FlipSelector_Check_v5
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
    Expected_Flip_Effect_v5 expFlipEffect
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
void Check_Given_Flip_Is_Optimal_v5
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


// Designed to be called at the very end of flips
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
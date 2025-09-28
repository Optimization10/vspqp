#ifndef FM_FLIPS_V7_H
#define FM_FLIPS_V7_H

#include "shared_api.h"
#include "ds_api.h"
#include "flips_common_api.h"

/*----------------------------------------------------------------------------*/
// Macros
/*----------------------------------------------------------------------------*/

#define IS_BETTER_POINT(accumGain, bestAccumGain, wgtImb, bestWgtImb) \
  ((accumGain) > (bestAccumGain) || \
  ((accumGain) == (bestAccumGain) && (wgtImb) <= (bestWgtImb)))

#define NO_GAIN_FLAG -gamma * n - gamma


/*----------------------------------------------------------------------------*/
// Types
/*----------------------------------------------------------------------------*/

typedef struct
{
    VSPINT var ; // 0 or 1
    VSPINT ind ; // 0....n-1
} Flip_Record ;


/*----------------------------------------------------------------------------*/
// Function Declarations
/*----------------------------------------------------------------------------*/
static void Select_Flip
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSP_problem          *problem,
    VSPINT            *retFlipVar,
    VSPINT            *retIndFlip,
    VSPINT           *retFlipGain
) ;

static void Do_Flip
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSPINT                    var, // variable of the flip
    VSPINT                      i, // index of the flip
    VSP_problem          *problem
) ;

static void Update_Local_Data_With_Flip
(
    VSPINT            *vars[2],
    VSPINT           *Avars[2],
    VSPINT            *wgts[2],
    VSPINT                   m,
    VSPINT             *nflips,
    VSPINT                  *H,
    VSPINT          *accumGain,
    VSPINT          *bestPoint,
    VSPINT      *bestAccumGain,
    VSPINT         *bestWgtImb,
    VSPINT             indFlip,
    VSPINT             varFlip,
    VSPINT            gainFlip,
    VSP_state         **CPs[2],
    Flip_Record        *memory,
    VSP_problem       *problem
) ;


// Initialize sized doubly linked lists and Ts. Ignore locked vertices. If
// randomize is true, then build the lists in a randomized manner by going
// through the vertices in a random order. Otherwise, go through the vertices
// in numerical order
// locked - If null, then ignored. One or more of its entries can also be NULL,
//          in which case that particular array is ignored
static void Initialize_Lists_And_Ts
(
    VSPINT                 *vars[2],
    VSPINT                *Avars[2],
    VSPINT               *locked[2],
    Doubly_LLs_Sized   *lists[2][2],
    VSPINT                 Ts[2][2],
    VSP_problem            *problem,
    VSPINT                randomize  // 0 or 1
) ;

static void Locked_Arr_Reset(VSPINT *locked, VSPINT n) ;

// Memory management:
// At the end, the data inside retStaes may be swapped with some checkpoint
// in 0, H, so it will be freed, or it will not be swapped in which case
// we don't want to free it and it won't be freed


// TODO: The following is a better rollback mechanism: Save a copy of the
// original point, and  keep an array xflipsBestPoint of size n.
// The ith entry is the number of flips done in xi between the original point
// until the current best point. Each time a new best point is found, we
// go back thorugh the memory and update this array.
// At the end, when we wish to roll back, we go to the copy of the original
// and apply flips where the number is odd.
// NOTE (IMPORTANT) In order to avoid saving a copy of the original,
// we can do this kind of backwards and maintain this kind of "diff" array
// with respect to the current point, not to the original. Not sure exactly
// how it would work but it's something to think about.
static void Rollback
(
    VSP_state            **CPs[2],
    VSP_state       *retStates[2],
    VSPINT                      m,
    VSPINT                      H,
    VSPINT              bestPoint,
    VSPINT                 nflips,
    Flip_Record           *memory,
    VSP_problem          *problem
) ;


static void Check_Data
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSP_problem          *problem
) ;

#endif
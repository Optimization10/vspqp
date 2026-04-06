#ifndef DS_API_H
#define DS_API_H

#include "VSPQP.h"
#include "shared_api.h"

/*------------------------------------------------------------------------------
 Contents:
    Direct Access List
    Dynamic Stack
    Linked List
    Max Queue v1
    Min Queue v1
    Max Queue v2
    Penalties
    Linked Lists
    Heaps
    Restart Decision
    Dynamic Heaps
    Miscellaneous
------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Direct Access List */
/*----------------------------------------------------------------------------*/
typedef struct
{
    VSPINT *list     ;
    VSPINT *lookup   ;
    VSPINT  size     ;
    VSPINT  max_size ;
} Direct_Access_List ;

Direct_Access_List *Construct_Direct_Access_List(VSPINT max_size) ;
void Destruct_Direct_Access_List(Direct_Access_List *list) ;

void Direct_Access_List_Insert
(
    Direct_Access_List      *list,
    VSPINT              new_entry
) ;

void Direct_Access_List_Delete_By_Idx
(
    Direct_Access_List       *list,
    VSPINT              idx_to_del
) ;

void Direct_Access_List_Delete_By_Value
(
    Direct_Access_List       *list,
    VSPINT              val_to_del
) ;

VSPINT Direct_Access_List_Get_Rand_Entry(Direct_Access_List *List) ;

VSPINT Direct_Access_List_Get_Size(Direct_Access_List *List) ;

void Direct_Access_List_Reset
(
    Direct_Access_List *list
) ;

void Direct_Access_List_Print
(
    Direct_Access_List *list,
    VSPINT                *c, // Assumed that list contains indices into c
    char                *msg  // Ignored if NULL
) ;

void Direct_Access_List_Check
(
    Direct_Access_List *List
) ;



/*----------------------------------------------------------------------------*/
/* Levels */
/*----------------------------------------------------------------------------*/

Levels* Levels_copy(Levels* levels_in) ;



void destruct_levels
(
    Levels* levels
) ;

bool Levels_is_member
(
    Problem *problem,
    Levels  *levels
) ;

void Levels_select
(
    VSPINT       level_num,
    Levels  *levels
) ;

void Levels_select_top
(
    VSPINT n_top_levels,
    Levels *levels
) ;

void Levels_select_bot
(
    VSPINT n_bot_levels,
    Levels *levels
) ;

/*----------------------------------------------------------------------------*/
// Doubly Linked Lists
/*----------------------------------------------------------------------------*/

/*--- Unsized ----------------------------------------------------------------*/
typedef struct
{
    VSPINT *entries;
    VSPINT *fwd;
    VSPINT *rev;
    VSPINT  n;
} Doubly_LLs_Unsized;


Doubly_LLs_Unsized *Doubly_LLs_Unsized_Init(VSPINT n);
void Doubly_LLs_Unsized_Destroy(Doubly_LLs_Unsized *lls);
void Doubly_LLs_Unsized_Insert
(
    Doubly_LLs_Unsized     *lls,
    VSPINT              whichll,
    VSPINT                  vtx
);

void Doubly_LLs_Unsized_Delete
(
    Doubly_LLs_Unsized       *lls,
    VSPINT                whichll,
    VSPINT                    vtx
);

void Doubly_LLs_Unsized_Check
(
    Doubly_LLs_Unsized    *lls0,
    Doubly_LLs_Unsized    *lls1,
    VSPINT                   *z,
    VSPINT                  *Aw,
    VSPINT              *locked
);

void Doubly_LLs_Unsized_Reset(Doubly_LLs_Unsized *lls);

void Doubly_LLs_Unsized_See_Entry
(
    Doubly_LLs_Unsized      *lls,
    VSPINT               whichll,
    VSPINT                  *ret
) ;

VSPINT *Doubly_LLs_Unsized_Spread
(
    Doubly_LLs_Unsized *lls
) ;

/*--- Sized ----------------------------------------------------------------*/
typedef struct
{
    VSPINT *entries ; // Size n + 1
    VSPINT *sizes ;
    VSPINT *fwd ;     // Size n
    VSPINT *rev ;     // Size n
    VSPINT n ;
} Doubly_LLs_Sized ;


// Sized lists
Doubly_LLs_Sized *Doubly_LLs_Sized_Init(VSPINT n);
void Doubly_LLs_Sized_Destroy(Doubly_LLs_Sized *lls);
void Doubly_LLs_Sized_Insert(Doubly_LLs_Sized *lls, VSPINT whichll, VSPINT vtx);
void Doubly_LLs_Sized_Delete(Doubly_LLs_Sized *lls, VSPINT whichll, VSPINT vtx);
void Doubly_LLs_Sized_Check(Doubly_LLs_Sized *lls0, Doubly_LLs_Sized *lls1, VSPINT *z, VSPINT *Aw, VSPINT *locked);
void Doubly_LLs_Sized_Reset(Doubly_LLs_Sized *lls) ;
void Doubly_LLs_Sized_See_Entry
(
    Doubly_LLs_Sized *lls, VSPINT whichll, VSPINT *ret
) ;
VSPINT Doubly_LLs_Sized_Get_Size(Doubly_LLs_Sized *lls, VSPINT whichll) ;

/*--- Simple lists, not in a struct ------------------------------------------*/
void Doubly_Linked_Lists_Insert(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT which_ll, VSPINT to_insert, VSPINT n);
void Doubly_Linked_Lists_Delete(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT which_ll, VSPINT to_delete, VSPINT n);
bool Doubly_Linked_Lists_Is_Member_Specific_List(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT which_ll, VSPINT member, VSPINT n);
bool Doubly_Linked_Lists_Is_Member_Any_List(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT member, VSPINT n);
void Doubly_Linked_Lists_Spread_One_List(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT whichll, VSPINT *spread, VSPINT n);
VSPINT *Doubly_Linked_Lists_Spread(VSPINT *entries, VSPINT *fwd, VSPINT *rev, VSPINT n);
void Doubly_Linked_Lists_Check(VSPINT *entries0, VSPINT *fwd0, VSPINT *rev0, VSPINT *entries1, VSPINT *fwd1, VSPINT *rev1, VSPINT *z, VSPINT *Aw, VSPINT *locked, VSPINT n);
void Doubly_Linked_Lists_See_Entry
(
    VSPINT *entries,
    VSPINT  whichll,
    VSPINT     *ret
) ;


/******************************************************************************
 * MAX HEAP V0
 ******************************************************************************/
void MaxHeap_Heapify_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  numNodes,
    VSPINT  *c
);

void MaxHeap_Bubble_Down_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  numNodes,
    VSPINT  *c,
    VSPINT   p
);

void MaxHeap_Insert_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toInsert
);

void MaxHeap_Delete_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toDel
);

void MaxHeap_See_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topVal
);

void MaxHeap_Get_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topVal
);

void MaxHeap_See_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topKey
);

void MaxHeap_Get_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topKey
);

void MaxHeap_Check_v0
(
    VSPINT      *heap,
    VSPINT    *lookup,
    VSPINT   numNodes,
    VSPINT         *c
) ;

/******************************************************************************
 * MIN HEAP V0
 ******************************************************************************/
void MinHeap_Bubble_Down_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  numNodes,
    VSPINT  *c,
    VSPINT   p
);

void MinHeap_Heapify_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  numNodes,
    VSPINT  *c
);

void MinHeap_Insert_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toInsert
);

void MinHeap_Delete_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT   toDel
);

void MinHeap_See_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topVal
);

void MinHeap_Get_Top_Val_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topVal
);

void MinHeap_See_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *c,
    VSPINT   numNodes,
    VSPINT  *topKey
);

void MinHeap_Get_Top_Key_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT  *numNodes,
    VSPINT  *c,
    VSPINT  *topKey
);

void MinHeap_Check_v0
(
    VSPINT  *heap,
    VSPINT  *lookup,
    VSPINT   numNodes,
    VSPINT  *c
) ;



/*----------------------------------------------------------------------------*/
/* Heaps V3 */
/*----------------------------------------------------------------------------*/


void min_heap_add_v3 /* returns the new size of the heap */
(
    VSPINT *heap, /* size n + 1. containing indices into X */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT new, /* the new new */
    VSPINT *nheap /* number of elements in heap not counting new one */
) ;

void min_heap_delete_v3 /* return new size of heap */
(
    VSPINT *heap, /* indices into X, 1..n on input.  size nheap+1 */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT *nheap /* number of items in heap */
) ;

void min_heap_build_v1 (
    VSPINT*    heap, /* on input, an unsorted set of indices */
    VSPFLOAT*     x, /* numbers to sort */
    VSPINT        n  /* number of elements to build into the heap */
) ;

void partial_max_sort
(
    VSPINT      *w,/* indices of numbers to sort. Will sort in place */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT       theta, /* find the theta largest numbers */
    VSPINT       n  /* length of w */
) ;

void min_heap_build_v3
(
    VSPINT*    heap, /* on input, an unsorted set of indices */
    Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT        n  /* number of elements to build into the heap */
) ;



void destruct_restart_decision(restart_decision *dec) ;

/* If using the rule, and if it returns true, then it also increments
   dec->restarts_done. If using levels, objIn and objOut should be EMPTY */
bool restart_decide
(
    restart_decision       *dec,
    VSPFLOAT                  objIn, /* If not using rule, EMPTY */
    VSPFLOAT                 objOut, /* If not using rule, EMPTY */
    Problem                *problem  /* If using levels, this is ignored   */
) ;



/*----------------------------------------------------------------------------*/
/* Dynamic Heaps */
/*----------------------------------------------------------------------------*/
typedef struct
{
    VSPINT *heap ;
    VSPINT *ptrs ; // WARNING: MAY NOT BELOG TO THIS HEAP ONLY. MAY BE SHARED BETWEEN HEAPS
    VSPINT *c ;    // the keys
    VSPINT  numVals ;
    VSPINT  maxNumVals ;
    VSPINT     minOrMax ; // +1 for max heap, -1 for min heap
} Dynamic_Heap ;

Dynamic_Heap *Dynamic_Heap_Init
(
    VSPINT *c, VSPINT *ptrs, VSPINT maxNumVals, VSPINT minOrMax
) ;
void Dynamic_Heap_Destroy(Dynamic_Heap *Heap) ;
void Dynamic_Heap_Insert(Dynamic_Heap *Heap, VSPINT toInsert) ;
void Dynamic_Heap_Get_Top_Val(Dynamic_Heap *Heap, VSPINT *topVal) ;
void Dynamic_Heap_See_Top_Val(Dynamic_Heap *Heap, VSPINT *topVal) ;
void Dynamic_Heap_Delete(Dynamic_Heap *Heap, VSPINT toDel) ;
void Dynamic_Heap_Check(Dynamic_Heap *Heap) ;


/*----------------------------------------------------------------------------*/
// Heaps_v4
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
    VSPINT  n ; // Number of vertices at current level
}  Min_Heap_v4 ;

// Heap containing unique values of c
// Expects entries in 1, 2, ..., gamma
typedef struct
{   
    VSPINT *heap ;
    VSPINT *heapPtrs ;
    VSPINT  sizeHeap ;
    VSPINT  gamma ; // Max size of the heap
    VSPINT *c ;
    VSPINT  n ; // Number of vertices at current level
} Max_Heap_v4 ;

// Min Heap functions
Min_Heap_v4 *Min_Heap_Init_v4
(
    VSPINT gamma,
    VSPINT    *c,
    VSPINT     n  // Number of vertices at current level
) ;
void Min_Heap_Destroy_v4(Min_Heap_v4 *heap);
void Min_Heap_Check_v4(Min_Heap_v4 *Heap);
void Min_Heap_Bubble_Down_v4(Min_Heap_v4 *Heap, VSPINT p);
void Min_Heap_Bubble_Up_v4(Min_Heap_v4 *Heap, VSPINT p, VSPINT valToBubble);
void Min_Heap_Build_v4(Min_Heap_v4 *Heap);
void Min_Heap_Insert_v4(Min_Heap_v4 *Heap, VSPINT valToInsert);
void Min_Heap_Delete_v4(Min_Heap_v4 *Heap, VSPINT toDel);
VSPINT Min_Heap_See_Top_Key_v4(Min_Heap_v4 *Heap) ;

// Max Heap functions
Max_Heap_v4 *Max_Heap_Init_v4
(
    VSPINT gamma,
    VSPINT    *c,
    VSPINT     n  // number of vertices at current level
) ;
void Max_Heap_Destroy_v4(Max_Heap_v4 *heap);
void Max_Heap_Check_v4(Max_Heap_v4 *Heap);
void Max_Heap_Bubble_Down_v4(Max_Heap_v4 *Heap, VSPINT p);
void Max_Heap_Bubble_Up_v4(Max_Heap_v4 *Heap, VSPINT p, VSPINT valToBubble);
void Max_Heap_Build_v4(Max_Heap_v4 *Heap);
void Max_Heap_Insert_v4(Max_Heap_v4 *Heap, VSPINT valToInsert);
void Max_Heap_Delete_v4(Max_Heap_v4 *Heap, VSPINT toDel);
VSPINT Max_Heap_See_Top_Key_v4(Max_Heap_v4 *Heap) ;

/*----------------------------------------------------------------------------*/
// IdxGroupsByC
/*----------------------------------------------------------------------------*/
typedef enum
{
    MIN_HEAP_FLAG,
    MAX_HEAP_FLAG
} MinMax_Heap_Flag ;

typedef struct
{
    Min_Heap_v4 *minHeap;
    Max_Heap_v4 *maxHeap;
    VSPINT      **groups ;
    VSPINT      *sizes ;
    VSPINT      *maxSizes ; // kth entry is num indices in c with ci = k
    VSPINT      *ptrs ;     // -1 for entries which are not here
    VSPINT       gamma ;

    const VSPINT *c ;
    const VSPINT *z ;
    const VSPINT *Aw ;
    VSPINT targetzi ;
    VSPINT targetAwi ;
    char   var ;

    MinMax_Heap_Flag flag ; // -1 for min, 1 for max
} IdxGroupsByC ;

void IdxGroupsByC_See_Random_Idx_Given_C
(
    IdxGroupsByC    *Groups,
    VSPINT          randInt, // Used for randInt % sizeGroup
    VSPINT        givenCval,
    VSPINT       *idxReturn  // Return
) ;

IdxGroupsByC *IdxGroupsByC_Init
(
    VSPINT  gamma,
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT  nvtx,
    VSPINT  targetzi,
    VSPINT  targetAwi,
    VSPINT *z,
    VSPINT *Aw,
    char    var,
    MinMax_Heap_Flag flag
) ;

void IdxGroupsByC_Destroy
(
    IdxGroupsByC *Groups
) ;

void IdxGroupsByC_Insert
(
    IdxGroupsByC *Groups,
    VSPINT        idxToInsert
) ;

void IdxGroupsByC_Delete
(
    IdxGroupsByC *Groups,
    VSPINT        idxToDel
) ;

void IdxGroupsByC_See_Top_Cval_If_Feasible
(
    IdxGroupsByC *Groups,
    VSPINT        wgt,
    VSPINT        u,
    VSPINT       *cvalRet
) ;

void IdxGroupsByC_Get_Random_Idx_Top_C
(
    IdxGroupsByC *Groups,
    VSPINT        randInt,
    VSPINT       *idxReturn
) ;

void IdxGroupsByC_Num_Idx_With_Given_C
(
    IdxGroupsByC    *Groups,
    VSPINT        givenCval,
    VSPINT       *retNumIdx
) ;

void IdxGroupsByC_See_Random_Idx_Top_C
(
    IdxGroupsByC *Groups,
    VSPINT        randInt,
    VSPINT       *idxReturn
) ;

void IdxGroupsByC_Num_Idx_With_Top_C
(
    IdxGroupsByC *Groups,
    VSPINT       *retNumIdx
) ;

void IdxGroupsByC_Check
(
    IdxGroupsByC *Groups,
    VSPINT             n  // Size of z and Aw arrays
) ;


void IdxGroupsByC_See_Best_Feasible_Cval
(
    IdxGroupsByC    *Groups,
    VSPINT              wgt,
    VSPINT                u,
    bool        deep_search, // If false, only look at top. If true, search heap
    VSPINT         *cvalRet, // EMPTY if not feasible
    VSPINT               *z, // Used for debugging only
    VSPINT              *Aw, // Used for debugging only
    VSPINT                n, // Used for debugging only
    VSPINT               *c  // Used for debugging only
) ;

/*----------------------------------------------------------------------------*/
/* Miscellaneous */
/*----------------------------------------------------------------------------*/






#endif
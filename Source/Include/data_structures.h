#ifndef DS_API_H
#define DS_API_H
#include "core.h"
#include "VSPQP.h"


/*------------------------------------------------------------------------------
 Contents:
    Direct Access List
    Heaps V4 (Min/Max, unique c-values)
    IdxGroupsByC
------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/* Direct Access List                                                         */
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
    VSPINT                *c,
    char                *msg
) ;

void Direct_Access_List_Check
(
    Direct_Access_List *List
) ;


/*----------------------------------------------------------------------------*/
/* Max Heap V0                                                                */
/*----------------------------------------------------------------------------*/
void MaxHeap_Heapify_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c
) ;

void MaxHeap_Bubble_Down_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c,
    VSPINT    p
) ;

void MaxHeap_Insert_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT    toInsert
) ;

void MaxHeap_Delete_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT    toDel
) ;

void MaxHeap_See_Top_Val_v0
(
    VSPINT   *heap,
    VSPINT   *c,
    VSPINT    numNodes,
    VSPINT   *topVal
) ;

void MaxHeap_Get_Top_Val_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT   *topVal
) ;

void MaxHeap_See_Top_Key_v0
(
    VSPINT   *heap,
    VSPINT   *c,
    VSPINT    numNodes,
    VSPINT   *topKey
) ;

void MaxHeap_Get_Top_Key_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT   *topKey
) ;

void MaxHeap_Check_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c
) ;


/*----------------------------------------------------------------------------*/
/* Min Heap V0                                                                */
/*----------------------------------------------------------------------------*/
void MinHeap_Bubble_Down_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c,
    VSPINT    p
) ;

void MinHeap_Heapify_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c
) ;

void MinHeap_Insert_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT    toInsert
) ;

void MinHeap_Delete_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT    toDel
) ;

void MinHeap_See_Top_Val_v0
(
    VSPINT   *heap,
    VSPINT   *c,
    VSPINT    numNodes,
    VSPINT   *topVal
) ;

void MinHeap_Get_Top_Val_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT   *topVal
) ;

void MinHeap_See_Top_Key_v0
(
    VSPINT   *heap,
    VSPINT   *c,
    VSPINT    numNodes,
    VSPINT   *topKey
) ;

void MinHeap_Get_Top_Key_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT   *numNodes,
    VSPINT   *c,
    VSPINT   *topKey
) ;

void MinHeap_Check_v0
(
    VSPINT   *heap,
    VSPINT   *lookup,
    VSPINT    numNodes,
    VSPINT   *c
) ;

/*----------------------------------------------------------------------------*/
/* Heaps v1 (same as v0, except no lookup)                                    */
/*----------------------------------------------------------------------------*/

void MaxHeap_Heapify_v1
(
    VSPINT       *heap,
    VSPINT    numNodes,
    VSPINT          *c
) ;

void MaxHeap_Insert_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT  toInsert
) ;

void MaxHeap_Delete_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   address
) ;

void MaxHeap_See_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topVal
) ;

void MaxHeap_Get_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topVal
) ;

void MaxHeap_See_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topKey
) ;

void MaxHeap_Get_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topKey
) ;

void MaxHeap_Check_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
) ;

/*----------------------------------------------------------------------------*/
/* Min heap                                                                    */
/*----------------------------------------------------------------------------*/

void MinHeap_Heapify_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
) ;

void MinHeap_Insert_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT  toInsert
) ;

void MinHeap_Delete_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   address
) ;

void MinHeap_See_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topVal
) ;

void MinHeap_Get_Top_Val_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topVal
) ;

void MinHeap_See_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT        *c,
    VSPINT  numNodes,
    VSPINT   *topKey
) ;

void MinHeap_Get_Top_Key_v1
(
    VSPINT     *heap,
    VSPINT *numNodes,
    VSPINT        *c,
    VSPINT   *topKey
) ;

void MinHeap_Check_v1
(
    VSPINT     *heap,
    VSPINT  numNodes,
    VSPINT        *c
) ;

/*----------------------------------------------------------------------------*/
// Heap v2 (with lookup) prototypes
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Heap v2 structs and prototypes
/*----------------------------------------------------------------------------*/

typedef struct
{
    VSPINT     *heap ;   /* 1-indexed, slots 1..n (owned) */
    VSPINT   *lookup ;   /* lookup[i] = heap position of i, or EMPTY (owned) */
    VSPINT        *c ;   /* key array (borrowed) */
    VSPINT         n ;   /* length of c = capacity */
    VSPINT  numNodes ;   /* current size */
} max_heap_v2 ;

typedef struct
{
    VSPINT     *heap ;
    VSPINT   *lookup ;
    VSPINT        *c ;
    VSPINT         n ;
    VSPINT  numNodes ;
} min_heap_v2 ;

/* Max heap */
VSPINT MaxHeap_Init_v2        (max_heap_v2 *h, VSPINT *c, VSPINT n) ;
void   MaxHeap_Free_v2        (max_heap_v2 *h) ;
void   MaxHeap_Reset_v2       (max_heap_v2 *h) ;
VSPINT MaxHeap_Contains_v2    (const max_heap_v2 *h, VSPINT i) ;
void   MaxHeap_Heapify_v2     (max_heap_v2 *h) ;
void   MaxHeap_Build_v2       (max_heap_v2 *h, const VSPINT *elems, VSPINT numElems) ;
void   MaxHeap_Insert_v2      (max_heap_v2 *h, VSPINT toInsert) ;
void   MaxHeap_Delete_v2      (max_heap_v2 *h, VSPINT address) ;
void   MaxHeap_Delete_Val_v2  (max_heap_v2 *h, VSPINT toDel) ;
void   MaxHeap_See_Top_Val_v2 (const max_heap_v2 *h, VSPINT *topVal) ;
void   MaxHeap_Get_Top_Val_v2 (max_heap_v2 *h, VSPINT *topVal) ;
void   MaxHeap_See_Top_Key_v2 (const max_heap_v2 *h, VSPINT *topKey) ;
void   MaxHeap_Get_Top_Key_v2 (max_heap_v2 *h, VSPINT *topKey) ;
void   MaxHeap_Check_v2       (const max_heap_v2 *h) ;

/* Min heap */
VSPINT MinHeap_Init_v2        (min_heap_v2 *h, VSPINT *c, VSPINT n) ;
void   MinHeap_Free_v2        (min_heap_v2 *h) ;
void   MinHeap_Reset_v2       (min_heap_v2 *h) ;
VSPINT MinHeap_Contains_v2    (const min_heap_v2 *h, VSPINT i) ;
void   MinHeap_Heapify_v2     (min_heap_v2 *h) ;
void   MinHeap_Build_v2       (min_heap_v2 *h, const VSPINT *elems, VSPINT numElems) ;
void   MinHeap_Insert_v2      (min_heap_v2 *h, VSPINT toInsert) ;
void   MinHeap_Delete_v2      (min_heap_v2 *h, VSPINT address) ;
void   MinHeap_Delete_Val_v2  (min_heap_v2 *h, VSPINT toDel) ;
void   MinHeap_See_Top_Val_v2 (const min_heap_v2 *h, VSPINT *topVal) ;
void   MinHeap_Get_Top_Val_v2 (min_heap_v2 *h, VSPINT *topVal) ;
void   MinHeap_See_Top_Key_v2 (const min_heap_v2 *h, VSPINT *topKey) ;
void   MinHeap_Get_Top_Key_v2 (min_heap_v2 *h, VSPINT *topKey) ;
void   MinHeap_Check_v2       (const min_heap_v2 *h) ;

/*----------------------------------------------------------------------------*/
/* Heaps V3 (min heap, ratio-keyed Az_i / c_i)                               */
/*----------------------------------------------------------------------------*/
void min_heapify_v3
(
    VSPINT   *heap,
    Az       *Az,
    VSPINT   *c,
    VSPINT    p,
    VSPINT    nheap
) ;

void min_heap_build_v3
(
    VSPINT   *heap,
    Az       *Az,
    VSPINT   *c,
    VSPINT    n
) ;

void min_heap_add_v3
(
    VSPINT   *heap,
    Az       *Az,
    VSPINT   *c,
    VSPINT    new,
    VSPINT   *nheap
) ;

void min_heap_delete_v3
(
    VSPINT   *heap,
    Az       *Az,
    VSPINT   *c,
    VSPINT   *nheap
) ;

/*----------------------------------------------------------------------------*/
/* Heaps V4 (unique c-values, min and max)                                    */
/*----------------------------------------------------------------------------*/
typedef struct
{
    VSPINT *heap ;
    VSPINT *heapPtrs ;
    VSPINT  sizeHeap ;
    VSPINT  gamma ;
    VSPINT *c ;
    VSPINT  n ;
} Min_Heap_v4 ;

typedef struct
{
    VSPINT *heap ;
    VSPINT *heapPtrs ;
    VSPINT  sizeHeap ;
    VSPINT  gamma ;
    VSPINT *c ;
    VSPINT  n ;
} Max_Heap_v4 ;

Min_Heap_v4 *Min_Heap_Init_v4(VSPINT gamma, VSPINT *c, VSPINT n) ;
void Min_Heap_Destroy_v4(Min_Heap_v4 *heap) ;
void Min_Heap_Check_v4(Min_Heap_v4 *Heap) ;
void Min_Heap_Bubble_Down_v4(Min_Heap_v4 *Heap, VSPINT p) ;
void Min_Heap_Bubble_Up_v4(Min_Heap_v4 *Heap, VSPINT p, VSPINT valToBubble) ;
void Min_Heap_Build_v4(Min_Heap_v4 *Heap) ;
void Min_Heap_Insert_v4(Min_Heap_v4 *Heap, VSPINT valToInsert) ;
void Min_Heap_Delete_v4(Min_Heap_v4 *Heap, VSPINT toDel) ;
VSPINT Min_Heap_See_Top_Key_v4(Min_Heap_v4 *Heap) ;

Max_Heap_v4 *Max_Heap_Init_v4(VSPINT gamma, VSPINT *c, VSPINT n) ;
void Max_Heap_Destroy_v4(Max_Heap_v4 *heap) ;
void Max_Heap_Check_v4(Max_Heap_v4 *Heap) ;
void Max_Heap_Bubble_Down_v4(Max_Heap_v4 *Heap, VSPINT p) ;
void Max_Heap_Bubble_Up_v4(Max_Heap_v4 *Heap, VSPINT p, VSPINT valToBubble) ;
void Max_Heap_Build_v4(Max_Heap_v4 *Heap) ;
void Max_Heap_Insert_v4(Max_Heap_v4 *Heap, VSPINT valToInsert) ;
void Max_Heap_Delete_v4(Max_Heap_v4 *Heap, VSPINT toDel) ;
VSPINT Max_Heap_See_Top_Key_v4(Max_Heap_v4 *Heap) ;


/*----------------------------------------------------------------------------*/
/* IdxGroupsByC                                                               */
/*----------------------------------------------------------------------------*/
typedef enum
{
    MIN_HEAP_FLAG,
    MAX_HEAP_FLAG
} MinMax_Heap_Flag ;

typedef struct
{
    Min_Heap_v4 *minHeap ;
    Max_Heap_v4 *maxHeap ;
    VSPINT      **groups ;
    VSPINT       *sizes ;
    VSPINT       *maxSizes ;
    VSPINT       *ptrs ;
    VSPINT        gamma ;

    const VSPINT *c ;
    const VSPINT *z ;
    const VSPINT *Aw ;
    VSPINT targetzi ;
    VSPINT targetAwi ;
    char   var ;

    MinMax_Heap_Flag flag ;
} IdxGroupsByC ;

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

void IdxGroupsByC_Destroy(IdxGroupsByC *Groups) ;

void IdxGroupsByC_Insert(IdxGroupsByC *Groups, VSPINT idxToInsert) ;

void IdxGroupsByC_Delete(IdxGroupsByC *Groups, VSPINT idxToDel) ;

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

void IdxGroupsByC_See_Random_Idx_Top_C
(
    IdxGroupsByC *Groups,
    VSPINT        randInt,
    VSPINT       *idxReturn
) ;

void IdxGroupsByC_See_Random_Idx_Given_C
(
    IdxGroupsByC *Groups,
    VSPINT        randInt,
    VSPINT        givenCval,
    VSPINT       *idxReturn
) ;

void IdxGroupsByC_Num_Idx_With_Top_C
(
    IdxGroupsByC *Groups,
    VSPINT       *retNumIdx
) ;

void IdxGroupsByC_Num_Idx_With_Given_C
(
    IdxGroupsByC *Groups,
    VSPINT        givenCval,
    VSPINT       *retNumIdx
) ;

void IdxGroupsByC_Check(IdxGroupsByC *Groups, VSPINT n) ;

void IdxGroupsByC_See_Best_Feasible_Cval
(
    IdxGroupsByC *Groups,
    VSPINT        wgt,
    VSPINT        u,
    bool          deep_search,
    VSPINT       *cvalRet,
    VSPINT       *z,
    VSPINT       *Aw,
    VSPINT        n,
    VSPINT       *c
) ;


#endif
#ifndef DS_API_H
#define DS_API_H

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
/* Dynamic Stack */
/*----------------------------------------------------------------------------*/


VSP_dyn_stk *VSP_construct_dyn_stk() ;
void VSP_destruct_dyn_stk(VSP_dyn_stk* ds) ;
void VSP_dyn_stk_push(VSP_dyn_stk* ds, VSPINT to_push) ;
bool VSP_dyn_stk_is_member(VSPINT i, VSP_dyn_stk* ds) ;
VSP_dyn_stk *VSP_dyn_stk_copy(VSP_dyn_stk *stack_in) ;

/*----------------------------------------------------------------------------*/
/* Levels */
/*----------------------------------------------------------------------------*/

VSP_levels* VSP_levels_copy(VSP_levels* levels_in) ;

VSP_levels* VSP_construct_levels
(
    VSP_levels_selection_type selection_type
) ;

void VSP_destruct_levels
(
    VSP_levels* levels
) ;

bool VSP_levels_is_member
(
    VSP_problem *problem,
    VSP_levels  *levels
) ;

void VSP_levels_select
(
    VSPINT       level_num,
    VSP_levels  *levels
) ;

void VSP_levels_unselect
(
    VSPINT     level_num,
    VSP_levels   *levels
) ;

void VSP_levels_select_top
(
    VSPINT n_top_levels,
    VSP_levels *levels
) ;

void VSP_levels_select_bot
(
    VSPINT n_bot_levels,
    VSP_levels *levels
) ;


/*----------------------------------------------------------------------------*/
/* Penalties */
/*----------------------------------------------------------------------------*/

typedef struct
{
    VSPINT *z0_entry ;    /* Entry points into the z0 linked lists */
    VSPINT *z0_fwd_list ;
    VSPINT *z0_rev_list ;

    VSPINT *z1_entry ;     /* Entry points into the z1 linked lists */
    VSPINT *z1_fwd_list ;
    VSPINT *z1_rev_list ;

    VSPINT max_pen ; /* max {Awi : zi = 1 }*/
    VSPINT min_pen ; /* min {Awi : zi = 0 }*/

    VSPINT n ;             /* Size of the matrix */

} VSP_penalties ;

VSP_penalties *VSP_construct_penalties(VSPINT n) ;
void VSP_destruct_penalties(VSP_penalties *pens) ;

void VSP_copy_penalties
(
    VSP_penalties *from,
    VSP_penalties *to
) ;

void VSP_build_penalties
(
    VSPINT                 *z,
    VSP_penalties     *Z_pens,
    VSPINT                *Aw,
    VSP_problem      *problem
) ;


void VSP_check_penalties
(
    VSP_state *X,
    VSP_state *Y,
    VSP_penalties *X_pen,
    VSP_penalties *Y_pen,
    VSPINT n
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
/* Max Queue v1 */
/*----------------------------------------------------------------------------*/
/* "priorities" is a 0-based array of size n. priorities[0 ... n-1]
   "queue" holds indices into the priorities array in queue[1 .. n_queue]
   The priority of i is data[i]. The highest data value has the highest
   priority (is at the root)                                                  */
typedef struct
{
    VSPINT *queue ;
    VSPINT *priorities ;
    VSPINT  n_queue ;      /* number of items currently in the queue */
    VSPINT  n_priorities ; /* number of items in the priorities array. Also,
                              this is the max number of items in the queue */
    VSPINT *lookup ;       /* lookup[i] = k means queue[k] = i. if k is EMPTY,
                              i is not in the queue */
} Max_Queue_v1 ;


Max_Queue_v1 *VSP_construct_max_queue_v1
(
    VSPINT n_priorities,
    VSPINT *priorities
) ;

void VSP_destruct_max_queue_v1
(
    Max_Queue_v1 *to_destroy
) ;

void VSP_max_queue_build_v1
(
    Max_Queue_v1 *queuestr
) ;

VSPINT Max_Queue_See_Top_Key_v1
(
    Max_Queue_v1 *queue
) ;

VSPINT VSP_max_queue_see_top_priority_v1
(
    Max_Queue_v1 *queuestr
) ;

void Max_Queue_Delete_v1
(
    Max_Queue_v1 *queuestr,
    VSPINT              rm
) ;

void Max_Queue_Reset_v1
(
    Max_Queue_v1 *qstr
) ;

void Max_Queue_Insert_v1
(
    Max_Queue_v1      *qstr,
    VSPINT        to_insert
) ;
/*----------------------------------------------------------------------------*/
/* Min Queue v1 */
/*----------------------------------------------------------------------------*/

/* "priorities" is a 0-based array of size n. priorities[0 ... n-1]
   "queue" holds indices into the priorities array in queue[1 .. n_queue]
   The priority of i is data[i]. The lowest data value has the highest
   priority (is at the root)                                                  */
typedef struct
{
    VSPINT *queue ;
    VSPINT *priorities ;
    VSPINT  n_queue ;      /* number of items currently in the queue */
    VSPINT  n_priorities ; /* number of items in the priorities array. Also,
                              this is the max number of items in the queue */
    VSPINT *lookup ;       /* lookup[i] = k means queue[k] = i. if k is EMPTY,
                              i is not in the queue */
} Min_Queue_v1 ;



VSPINT VSP_min_queue_see_top_entry_v1
(
    Min_Queue_v1 *queuestr
) ;

VSPINT VSP_min_queue_see_top_priority_v1
(
    Min_Queue_v1 *queuestr
) ;


/* toadd is an index into the priorities array */
void VSP_min_queue_add_v1
(
    Min_Queue_v1 *queuestr,
    VSPINT toadd
) ;

/* Remove entry rm in the queue. rm is an index of the priorities array, not
   an index of the heap array. */
void Min_Queue_Delete_v1
(
    Min_Queue_v1 *queuestr,
    VSPINT rm
) ;

/* This function heapifies the queue array. It assumes the lookup array is 
   correct */
void VSP_min_queue_build_v1
(
    Min_Queue_v1 *queuestr
) ;

/* Construct an uninitialized min priority queue v1 ds. Everything including the
   lookup is initialized */
Min_Queue_v1 *VSP_construct_min_queue_v1
(
    VSPINT n_priorities,
    VSPINT *priorities
) ;

void VSP_destruct_min_queue_v1
(
    Min_Queue_v1 *to_destroy
) ;

void Min_Queue_Reset_v1
(
    Min_Queue_v1 *qstr
) ;

// ASSUMED: to_insert is NOT in the queue
void Min_Queue_Insert_v1
(
    Min_Queue_v1      *qstr,
    VSPINT        to_insert
) ;

void Min_Queue_Print_Node_v1
(
    Min_Queue_v1   *q,
    VSPINT       node,
    char         *msg  // Ignored if NULL
) ;

void Max_Queue_Print_Node_v1
(
    Max_Queue_v1   *q,
    VSPINT       node,
    char         *msg  // Ignored if NULL
) ;

void VSP_check_min_heap_v1
(
    Min_Queue_v1* queuestr
) ;


/*----------------------------------------------------------------------------*/
/* Max Queue v2 */
/*----------------------------------------------------------------------------*/

/* "priorities" is a 0-based array of size n. priorities[0 ... n-1]
   "queue" holds indices into the priorities array in queue[1 .. n_queue]
   The priority of i is top[i] / bot[i]. The largest fraction has the highest
   priority (is at the root)                                                  */
typedef struct
{
    VSPINT *queue ;
    VSPINT  *top ;         /* numerators of priorities */
    VSPINT  *bot ;         /* denominators of priorities */
    VSPINT  n_queue ;      /* number of items currently in the queue */
    VSPINT  n_priorities ; /* number of items in the priorities array. Also,
                              this is the max number of items in the queue */
    VSPINT *lookup ;       /* lookup[i] = k means queue[k] = i. if k is EMPTY,
                              i is not in the queue */
} VSP_max_queue_v2 ;


VSPFLOAT VSP_max_queue_see_top_priority_v2
(
    VSP_max_queue_v2 *queuestr
) ;

VSP_max_queue_v2 *VSP_construct_max_queue_v2
(
    VSPINT n_priorities, /* max size of the queue */
    VSPINT *top,
    VSPINT *bot
) ;

void VSP_max_queue_add_v2
(
    VSP_max_queue_v2 *queuestr,
    VSPINT toadd
) ;

void VSP_destruct_max_queue_v2
(
    VSP_max_queue_v2 *to_destroy
) ;

void VSP_max_queue_build_v2
(
    VSP_max_queue_v2 *queuestr
) ;

void VSP_max_queue_delete_v2
(
    VSP_max_queue_v2 *queuestr,
    VSPINT rm
) ;


VSPINT VSP_max_queue_see_top_entry_v2
(
    VSP_max_queue_v2 *queuestr
) ;

VSPFLOAT VSP_max_queue_see_top_priority_v2
(
    VSP_max_queue_v2 *queuestr
) ;

void VSP_max_queue_delete_wo_check_v2
(
    VSP_max_queue_v2 *queuestr,
    VSPINT rm
) ;

void VSP_max_queue_reset_v2
(
    VSP_max_queue_v2 *q
) ;



bool VSP_max_queue_v2_is_member
(
    VSPINT i,
    VSP_max_queue_v2 *q
) ;

void VSP_max_queue_v2_bubble_up
(
    VSP_max_queue_v2 *q,
    VSPINT node
) ;


/*----------------------------------------------------------------------------*/
/* Heaps */
/*----------------------------------------------------------------------------*/


void VSP_min_heap_add_v2 /* returns the new size of the heap */
(
    VSPINT *heap, /* size n + 1. containing indices into X */
    VSPINT *num,  /* size n. Numerator vector to sort by */
    VSPINT *denom,  /* size n. Denominator vector to sort by */
    VSPINT new, /* the new new */
    VSPINT *nheap /* number of elements in heap not counting new one */
) ;

void VSP_check_max_heap_v2
(
    VSP_max_queue_v2* queuestr,
    VSPFLOAT               gamma
) ;

void VSP_max_heapify_v2
(
    VSP_max_queue_v2 *queuestr,
    VSPINT            p          /* start at node p in the heap */
) ;

void VSP_min_heap_delete_v2 /* return new size of heap */
(
    VSPINT *heap, /* indices into X, 1..n on input.  size nheap+1 */
    VSPINT *num,  /* size n. Numerator vector to sort by */
    VSPINT *denom,  /* size n. Denominator vector to sort by */
    VSPINT *nheap /* number of items in heap */
) ;

void VSP_min_heap_build_v2
(
    VSPINT*    heap, /* on input, an unsorted set of indices */
    VSPINT*     num, /* size n. Numerator vector to sort by */
    VSPINT*     denom, /* size n. Denominator vector to sort by */
    VSPINT        n  /* number of elements to build into the heap */
) ;


void VSP_min_heap_add_v3 /* returns the new size of the heap */
(
    VSPINT *heap, /* size n + 1. containing indices into X */
    VSP_Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT new, /* the new new */
    VSPINT *nheap /* number of elements in heap not counting new one */
) ;

void VSP_min_heap_delete_v3 /* return new size of heap */
(
    VSPINT *heap, /* indices into X, 1..n on input.  size nheap+1 */
    VSP_Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT *nheap /* number of items in heap */
) ;

void VSP_min_heap_build_v1 (
    VSPINT*    heap, /* on input, an unsorted set of indices */
    VSPFLOAT*     x, /* numbers to sort */
    VSPINT        n  /* number of elements to build into the heap */
) ;

void VSP_partial_max_sort
(
    VSPINT      *w,/* indices of numbers to sort. Will sort in place */
    VSP_Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT       theta, /* find the theta largest numbers */
    VSPINT       n  /* length of w */
) ;

void VSP_min_heap_build_v3
(
    VSPINT*    heap, /* on input, an unsorted set of indices */
    VSP_Az*     Az, /* size n. Numerator vector to sort by */
    VSPINT*     c, /* size n. Denominator vector to sort by */
    VSPINT        n  /* number of elements to build into the heap */
) ;


/*----------------------------------------------------------------------------*/
/* Restart Decision */
/*----------------------------------------------------------------------------*/



VSP_restart_decision *VSP_construct_restart_decision() ;
void VSP_destruct_restart_decision(VSP_restart_decision *dec) ;

/* If using the rule, and if it returns true, then it also increments
   dec->restarts_done. If using levels, objIn and objOut should be EMPTY */
bool VSP_restart_decide
(
    VSP_restart_decision       *dec,
    VSPFLOAT                  objIn, /* If not using rule, EMPTY */
    VSPFLOAT                 objOut, /* If not using rule, EMPTY */
    VSP_problem            *problem  /* If using levels, this is ignored   */
) ;

void VSP_restart_decide_reset(VSP_restart_decision *dec) ;


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
/*----------------------------------------------------------------------------*/
/* Miscellaneous */
/*----------------------------------------------------------------------------*/






#endif
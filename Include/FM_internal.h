#ifndef FM_INTERNAL_H
#define FM_INTERNAL_H

#include "shared_api.h"
#include "FM_api.h" // IWYU pragma: keep


#define INC_DEC(a, b, val) do {(a) += (val); (b) -= (val);} while(0)

/*----------------------------------------------------------------------------*/
/*      Types                                                                 */
/*----------------------------------------------------------------------------*/
typedef struct
{
    VSPINT value ; // Nodes
    VSPINT key ;
} Node ;

typedef struct
{
    Node   *heap        ;
    VSPINT *locator     ;   
    VSPINT numNodes     ;
    VSPINT maxNumNodes  ;
} QueueFM ;

/* Each node in the sep is associated with one of these structs.
   wgts[0] is the sum of weights (costs) of the neighbors in set 0 (i.e. x)
   and wgts[1] is the sum of weigts of the neighbors in set 1 (i.e. y)        */
typedef struct
{
    VSPINT wgts[2] ;
} ExtNeighWgts ;

/*----------------------------------------------------------------------------*/
/*      Function Prototypes                                                   */
/*----------------------------------------------------------------------------*/

void BND_Check
(
    VSPINT        nbnd,
    VSPINT     *bndind,
    VSPINT     *bndptr,
    VSPINT      *where,
    VSPINT           n  // size of where array
) ;


void QueueFM_Update(QueueFM *q, VSPINT value, VSPINT newKey) ;

VSPINT QueueFM_Get_Top_Value(QueueFM *q) ;

void QueueFM_Delete(QueueFM *q, VSPINT valueToDelete) ;

QueueFM *QueueFM_Construct
(
    VSPINT maxNumNodes
) ;

void QueueFM_Destruct
(
    QueueFM *q
) ;

void QueueFM_Check_Heap_Light(QueueFM *q) ;

void QueueFM_Insert(QueueFM *q, VSPINT newValue, VSPINT newKey) ;

VSPINT QueueFM_See_Top_Value(QueueFM *q) ;

VSPINT QueueFM_See_Top_Key(QueueFM *q) ;


#endif
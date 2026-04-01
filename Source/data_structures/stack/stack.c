/*----------------------------------------------------------------------------*/
// Simple implementation of a stack, WITHOUT indexing.
/*----------------------------------------------------------------------------*/

#include "stack.h"

struct StackInt
{
    VSPINT *stack ;
    VSPINT  top ;
    VSPINT  capacity ;
} ;

StackInt *StackInt_Init(VSPINT capacity)
{
    StackInt *ret = malloc(sizeof(StackInt)) ;
    ret->capacity = capacity ;
    ret->top   = -1 ;
    ret->stack = malloc(capacity * sizeof(VSPINT)) ;
    return ret ;
}

void StackInt_Destroy(StackInt *Stack)
{
    free(Stack->stack) ;
    free(Stack) ;
}

VSPINT StackInt_Pop(StackInt *Stack)
{
    if(Stack->top == -1) return EMPTY ;
    else return Stack->stack[Stack->top--] ;
}

// ASSUMED: we do not go over capacity
void StackInt_Push(StackInt *Stack, VSPINT item)
{
    Stack->stack[++Stack->top] = item ;
    return ;
}

VSPINT StackInt_GetSize(StackInt *Stack)
{
    return Stack->top + 1 ;
}
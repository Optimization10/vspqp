#include "dynamic_stack.h"

struct DynamicStack
{
    void **stack ;
    VSPINT size ;
    VSPINT maxSize ;
} ;

DynamicStack *DynamicStack_Init()
{
    VSPINT maxNumObjects = 10 ;
    DynamicStack *ret = calloc(1, sizeof(DynamicStack) ) ;
    ret->stack = calloc(maxNumObjects, sizeof(void*) ) ;
    ret->size = 0 ;
    ret->maxSize = maxNumObjects ;
    return ret ;
}

void DynamicStack_Destroy(DynamicStack *Dstack)
{
    free(Dstack->stack) ;
    free(Dstack) ;
}

void *DynamicStack_Pop(DynamicStack *Dstack)
{
    if( Dstack->size  == 0 ) return NULL ;
    else return Dstack->stack[ --Dstack->size ] ;
}

void DynamicStack_Push(DynamicStack *Dstack, void* newObject)
{
    if(Dstack->size == Dstack->maxSize)
    {
        Dstack->maxSize *= 2 ;
        void **temp = realloc(Dstack->stack, Dstack->maxSize * sizeof(void*) );
        if(!temp) return ;
        Dstack->stack = temp ;
    }

    Dstack->stack[ Dstack->size++ ] = newObject ;
    return ;
}

VSPINT DynamicStack_GetSize(DynamicStack *Dstack)
{
    return Dstack->size ;
}
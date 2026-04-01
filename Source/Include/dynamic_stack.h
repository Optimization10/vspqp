#ifndef DYNAMIC_STACK_H
#define DYNAMIC_STACK_H

#include "shared_api.h"

typedef struct DynamicStack DynamicStack ;


DynamicStack *DynamicStack_Init() ;
void DynamicStack_Destroy(DynamicStack *Dstack) ;
void DynamicStack_Push(DynamicStack *Dstack, void* newObject) ;
void *DynamicStack_Pop(DynamicStack *Dstack) ;
VSPINT DynamicStack_GetSize(DynamicStack *Dstack) ;

#endif
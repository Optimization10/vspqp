#ifndef STACK_H
#define STACK_H

#include "shared_api.h"

typedef struct StackInt StackInt ;

StackInt *StackInt_Init(VSPINT capacity) ;
void StackInt_Destroy(StackInt *Stack) ;
VSPINT StackInt_Pop(StackInt *Stack) ;
void StackInt_Push(StackInt *Stack, VSPINT item) ;
VSPINT StackInt_GetSize(StackInt *Stack) ;

#endif
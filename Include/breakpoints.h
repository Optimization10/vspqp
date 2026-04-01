#ifndef FM_FLIPS_V8_BREAKPOINTS_H
#define FM_FLIPS_V8_BREAKPOINTS_H
#include "shared_api.h"

typedef struct Breakpoints Breakpoints ;

Breakpoints *Breakpoints_Init(VSPINT lo, VSPINT hi, VSPINT numBP) ;
void Breakpoints_Destroy(Breakpoints *BP) ;
VSPINT Breakpoints_Get_NumBP(Breakpoints *BP) ;
VSPINT Breakpoints_Get_Point(Breakpoints *BP, VSPINT i) ;
VSPINT Breakpoints_Get_Next_Bp(Breakpoints *BP) ;

#endif
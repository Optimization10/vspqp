#include "shared_api.h"
#include "breakpoints.h"

/*------------------------------------------------------------------------*/
// breakpoints module
// purpose: provide equally spaced breakpoints between two integers lo and hi

// bp[0] = lo ---- bp[1] -------- bp[2] -------- .... ------- bp[numBP] = hi
/*------------------------------------------------------------------------*/

struct Breakpoints
{
    VSPINT lo ;
    VSPINT hi ;
    VSPINT *bp ;   // bp[0] = lo and bp[numBp] = hi
    VSPINT numBp ; // number of breakpoints including hi but excluding lo

    VSPINT bpCounter ; 
} ;


Breakpoints *Breakpoints_Init(VSPINT lo, VSPINT hi, VSPINT numBP)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(lo <= hi) ;
    VSP_ASSERT(numBP >= 1) ;

    Breakpoints *BP = malloc( sizeof(Breakpoints) ) ;
    BP->bp = malloc(sizeof(VSPINT) * (numBP + 1) ) ;
    BP->numBp = numBP ;
    BP->lo = lo ;
    BP->hi = hi ;


    VSPINT d = (hi - lo) / CAST(numBP, VSPFLOAT) ;
    for(VSPINT i = 0 ; i < numBP ; i++)
    {
        BP->bp[i] = lo + d * i ; // rounds down
        VSP_ASSERT(BP->bp[i] <= hi) ;
    }
    BP->bp[numBP] = hi ;
    return BP ;
}

void Breakpoints_Destroy(Breakpoints *BP)
{
    if(!BP) return ;
    free(BP->bp) ;
    free(BP) ;
}

VSPINT Breakpoints_Get_NumBP(Breakpoints *BP)
{
    return BP->numBp ;
}

// Get breakpoint i where i in [0 ... numBP]
VSPINT Breakpoints_Get_Point(Breakpoints *BP, VSPINT i)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT( IN_RANGE(i, 0, BP->numBp + 1) ) ;
    return BP->bp[i] ;
}

VSPINT Breakpoints_Get_Next_Bp(Breakpoints *BP)
{
    if(BP->bpCounter < BP->numBp)
    {
        return BP->bp[ (BP->bpCounter++) ] ;  
    }
    else
    {
        return BP->bp[BP->numBp] ;
    }
}
#include "core.h"
#include "FlipsCtxMain.h"
#include "ds_api.h"

struct FlipsCtxMain
{
    VSPINT *vars[2] ;
    VSPINT *Avars[2] ;
    Doubly_LLs_Sized *Lists[2][2] ;
    VSPINT Ts[2][2] ;
} ;
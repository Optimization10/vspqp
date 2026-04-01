#include "core.h"
#include "greedy_flips_ctx.h"
#include "ds_api.h"

struct GreedyFlipsCtx
{
    VSPINT *vars [2] ;
    VSPINT *Avars[2] ;
    VSPINT  wgts [2] ;
    VSPINT *locked[2] ;

    VSPINT Ts[2][2] ;
    VSPINT n ;
    VSP_problem *problem ;

    bool randomize ; // randomize initialization of lists
} ;
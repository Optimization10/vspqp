#include "VS.h"
#define EMPTY -1

VSPFLOAT *vs_to_bisect
(
    VSPINT   *shore_x,
    VSPINT    ind_frac_x,
    VSPINT   *shore_y,
    VSPINT    ind_frac_y,
    VSPINT    n
)
{
    VSPINT i ;
    VSPFLOAT *x = malloc(n * sizeof(*x)) ;

    for(i = 0 ; i < n ; i++)
    {
        if     (shore_x[i]) x[i] = 0   ;
        else if(shore_y[i]) x[i] = 1   ;
        else                x[i] = 0.5 ;
    }
    if(ind_frac_x != EMPTY) x[ind_frac_x] = 1 ;
    if(ind_frac_y != EMPTY) x[ind_frac_y] = 1 ;

    return x ;
}
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

    if(ind_frac_x != EMPTY) x[ind_frac_x] = 0 ;
    if(ind_frac_y != EMPTY) x[ind_frac_y] = 1 ;

    // Alternative. If val_frac_x is almost 0, then ind_frac_x is nearly
    // in the separator, so x[ind_frac_x] should be nearly 0.5. If val_frac_x
    // is almost 1, then x[ind_frac_x] should be nearly 0. Similar remarks for y
    // if(ind_frac_x != EMPTY) x[ind_frac_x] = 0.5 - val_frac_x / 2 ; 
    // if(ind_frac_y != EMPTY) x[ind_frac_y] = 0.5 + val_frac_y / 2 ; 

    return x ;
}
#include "core.h"
#include "hubs.h"

// Return a list of vertices whose degree / n is > t
void Get_Hubs
(
    VSPINT      n,
    VSPLONG    *Ap,
    VSPFLOAT    t,
    VSPINT    **hubs,    // Return
    VSPINT     *num_hubs // Return
)
{
    VSPINT i ;
    VSPFLOAT n_reciprocal = 1 / (VSPFLOAT) n ;

    *hubs     = malloc(n * sizeof(**hubs)) ;
    *num_hubs = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        if(Ap[i] * n_reciprocal > t)
        {
            (*hubs)[ (*num_hubs)++ ] = i ;
        }
    }
}
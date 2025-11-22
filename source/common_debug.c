#include "shared_api.h"

// ASSUMED: x, y are almost binary points  
void Check_Feasibility
(
    VSPINT      *xInt, // integer part of x 
    VSPINT   indFracX, // ignored if EMPTY
    VSPFLOAT valFracX,
    VSPINT      *yInt, // integer part of y
    VSPINT   indFracY, // ignored if EMPTY
    VSPFLOAT valFracY,
    VSPINT         *c,
    VSPINT          u,
    VSPINT          n
)
{
    DBG_FN_BOILERPLATE ;
    VSP_ASSERT(n > 0) ;
    VSP_ASSERT(indFracX == EMPTY || IN_RANGE(indFracX, 0, n) ) ;
    VSP_ASSERT(indFracY == EMPTY || IN_RANGE(indFracY, 0, n) ) ;

    // Check that the integer parts are binary and compute the weights
    VSPINT wgtx = 0, wgty = 0, i ;
    for(i = 0 ; i < n ; i++)
    {
        if(xInt[i] != 0 && xInt[i] != 1)
        {
            ERR("x not binary") ;
            break ;
        }

        if(yInt[i] != 0 && yInt[i] != 1)
        {
            ERR("y not binary") ;
            break ;
        }
        wgtx += xInt[i] * c[i] ;
        wgty += yInt[i] * c[i] ;
    }

    VSPFLOAT wgtFracX, wgtFracY ;
    if(indFracX != EMPTY)
    {
        wgtFracX = c[indFracX] * valFracX ;
        if( !VSP_equal(wgtFracX, lround(wgtFracX) ) )
        {
            ERR("WARNING: Weight of frac entry is not an integer") ;
        }
        wgtx    += (VSPINT) lround( wgtFracX ) ;
    }
    if(indFracY != EMPTY)
    {
        wgtFracY = c[indFracY] * valFracY ;
        if( !VSP_equal(wgtFracY, lround(wgtFracY) ) )
        {
            ERR("WARNING: Weight of frac entry is not an integer") ;
        }
        wgty    += (VSPINT) lround( wgtFracY ) ;
    }

    // Check that the weights satisfy the balance constraints
    if(wgtx < 1)
    {
        ERR("x fails lower balance constraint") ;
    }

    if(wgty < 1)
    {
        ERR("y fails lower balance constraint") ;
    }

    if(wgtx > u)
    {
        ERR("x fails upper balance constraint") ;
    }

    if(wgty > u)
    {
        ERR("y fails upper balance constraint") ;
    }

    return ;
}


// Assuming z is an integer vector, check that Az is correct
void Check_Az_Int
(
    const VSPINT *const  z,
    const VSPINT *const Az,
    const VSPINT *const Ai,
    const VSPLONG *const Ap,
    const VSPINT         n
)
{
    DBG_FN_BOILERPLATE ;
    VSP_ASSERT(n > 0) ;
    VSPINT *AzDbg = malloc(n * sizeof(*AzDbg)) ;
    VSP_mat_vec_mult(z, AzDbg, Ai, Ap, n) ;
    VSP_check_arrays_are_equal(Az, AzDbg, n) ;
    free(AzDbg) ;
}

void Check_Weight
(
    const VSPINT    *const z,
    const VSPINT     indFrac, // Ignored if EMPTY
    const VSPFLOAT   valFrac, // Ignored if EMPTY
    const VSPINT    *const c,
    const VSPINT         wgt,
    const VSPINT           n
)
{
    DBG_FN_BOILERPLATE ;
    VSP_ASSERT( indFrac == EMPTY || IN_RANGE(indFrac, 0, n) ) ;

    VSPINT       wgtDbg = 0 ;

    // Compute weight of integer part
    for (VSPINT i = 0 ; i < n ; i++)
    {
        wgtDbg += c[i] * z[i] ;
    }

    // Compute weight of fractional part
    if(indFrac != EMPTY)
    {
        wgtDbg += (VSPINT) lround( c[indFrac] * valFrac ) ;
    }

    if( wgtDbg != wgt )
    {
        VSP_error("Incorrect weight", DBG_INFO) ;
    }
    return ;
}


void Print_Neighs
(
    VSPINT vtx, VSPINT *Ai, VSPLONG *Ap
)
{
    FILE *file = fopen("output/Log_Dbg.txt", "a") ;
    fprintf(file, "\nNeighbors of vertex %lld :\n", (long long) vtx) ;
    for(VSPLONG p = Ap[vtx] ; p < Ap[vtx + 1] ; p++)
    {
        fprintf(file, "%lld  ", (long long) Ai[p]) ;
    }
    fclose(file) ;
    return ;
}
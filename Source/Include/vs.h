#ifndef VS_H
#define VS_H
#include "core.h"

typedef struct VSparm      /* parameters */
{
    /* printing */
    int       PrintLevel ; /* Level 0  = no print, ... , Level 2 = max print*/
    int       PrintParms ; /* TRUE means print parameter values */
    int      PrintStatus ; /* TRUE means print status */

    int      compute_eig ; /* TRUE means use eigenvector to start bisection */
    unsigned    int seed ; /* seed to start srand */
    VSPFLOAT         eps ; /* machine epsilon for VSFLOAT */
    int      SmallMatrix ; /* full diagonalization of A when n <= SmallMatrix */
    VSPFLOAT basis_factor ; /* # Lanczos iter = basis_factor*log2(dimension) */
    VSPFLOAT  Lanczos_tol ; /* |u_j| <= Lanczos_tol*Amax => stop Lanc. */
    VSPFLOAT       qr_its ; /* total QR iterations <= qr_its * dimension */
    VSPFLOAT   qr_tol_rel ; /* stop QR when u_{k-1}^2 <= qr_tol_rel*|d_k| */
    VSPFLOAT   qr_tol_abs ; /* stop QR when |u_{k-1}| <= qr_abs_tol */
    VSPFLOAT     qr_lower ; /* controls rescaling */
    int              ngp ; /* number of gradient projection iterations */
    int  initial_project ; /* TRUE: project z onto feasible set for bisection
                              before the gradient projection loop;
                              FALSE: skip the initial projection */

    int try_both_covers ;
} VSparm ;

void vs_default
(
    VSparm *Parm /* pointer to a vs parameter structure */
) ;

int vs
(
    int       *x, /* x [i] = 0, 1, or 2 to put i in set 0, set 1, or separator
                     respectively,  malloc'd by user */
    VSPINT    *Ei, /* row numbers */
    VSPLONG   *Ep, /* column pointers */
    VSPINT    *Ex, /* edge weights,   NULL => all entries in Ex are 1 */
    VSPINT     *a, /* vertex weights, NULL => all entries in a are 1 */
    VSPINT      n, /* problem dimension */
    VSPFLOAT    m, /* balance constraint parameter: 0.5 * original # vertices
                     m = 0.5*sum(vw [i]) */
    VSparm *Parm  /* parameter structure */
) ;

VSPINT VSPboundary
(
    VSPFLOAT           *x, /* feasible point for (GB), output described above */
    VSPFLOAT           *a, /* vertex weights, NULL => all entries in a are 1 */
    VSPFLOAT           *g, /* linear vector in objective */
    VSPINT              n, /* problem dimension */
    int const  PrintLevel
) ;

void VSP_Grad_Proj
(
    VSPFLOAT  *x,
    VSPFLOAT  *gx,
    VSPFLOAT  *y,
    VSPFLOAT  *gy,
    VSPFLOAT  *a,
    VSPINT    *Ai,
    VSPLONG   *Ap,
    VSPINT     u,
    VSPINT     n,
    VSPINT     gamma
) ;

#endif
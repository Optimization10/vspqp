#ifndef VS_H
#define VS_H
#include "VSPQP.h"
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

    /* vertex cover of bipartite graph */
    // First, select a cover type:
    // MIN_COST_COVER -- find a min cost cover by solving a max flow problem
    // MIN_CARD_COVER -- find a min cardinality using a max matching
    //
    // APPROX_MIN_CARD_COVER -- use the pricing algorithm to find a cover with
    //                          cardinality <= 2*OPT where OPT is the min
    //                          cardinality among all vertex covers
    cover_type cover_type ;

    // if cover_type = MIN_CARD_COVER, this boolean option controls whether
    //                 we compute a cover starting from each side and select
    //                 the cover with lower cost
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
    VSPFLOAT          *gx, /* objective gradient in subspace */
    VSPFLOAT          *gy, /* complementary objective gradient */
    VSPINT              n, /* problem dimension */
    VSPINT            *Si, /* row indices in subspace */
    VSPLONG           *Sp, /* col pointer in subspace */
    VSPFLOAT        gamma, /* penalty used in gradient */
    int const      update, /* = 1 if gy updated whenever x is updated, else 0 */
    int const  full_space,
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
    VSPFLOAT   gamma,
    int        ngp
) ;

// Only use this when the subproblem is much smaller than original problem
void VSP_FMsubspace
(
    VSPFLOAT  *x,    // variable in the bilinear objective
    VSPINT     fx,    // index of fractional component of x if one exists
    VSPFLOAT  *y,    // variable in the bilinear objective
    VSPINT    fy,    // index of fractional component of y if one exists
    VSPFLOAT  *a,    // vertex costs
    VSPINT    *Ai,   // row indices for full space
    VSPLONG   *Ap,   // col pointer for full space
    VSPINT     u,    // upper bound for balance constraint
    VSPINT     n,    // number of vertices
    VSPFLOAT   gamma,// max vertex cost
    VSPFLOAT   noFM  /*, skip FMsubspace when nF >= noFM*n
    VSPFLOAT   wx,      VSFLOAT workspace of size 7n
    VSPINT     wi,      VSINT workspace of size 3n + [nnz(A) +n]
    VSPLONG    wL       VSLINT workspace of size n+1 */
) ;

void VSP_Grad
(
    VSPFLOAT             *gx,
    VSPFLOAT       *const Gx,
    const VSPFLOAT *const y,
    const VSPINT          n,
    const VSPINT   *const Si,
    const VSPLONG  *const Sp,
    const VSPFLOAT        gamma,
    const int             FullSpace // TRUE if n for Full Space (not Subspace)
) ;

void FM_subspace_grad_proj(State *X, State *Y, Problem *P, Options *opts) ;

// Make part 2 a separator by moving a min-weight vertex cover of the
// bipartite graph between parts 0 and 1 into part 2. Defined in vs.c
int vs_separator_min_weight
(
    int      *part,
    VSPINT   *Ai,
    VSPLONG  *Ap,
    VSPINT    n,
    VSPINT   *vtx_wgts   // must not be NULL
) ;


void Grad_Proj_Full_Space
(
    State   *X,
    State   *Y,
    Problem *problem,
    Options *options
) ;

#endif
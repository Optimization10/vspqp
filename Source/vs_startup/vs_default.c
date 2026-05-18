#include "VS.h"

/* ==========================================================================
   ===== vs_default =========================================================
   ==========================================================================
    Set default parameter values.
   ========================================================================== */
void vs_default
(
    VSparm *Parm /* pointer to a vs parameter structure */
)
{
    /* printing */
    Parm->PrintLevel  = 0 ;
    Parm->PrintParms  = VSFALSE ;
    Parm->PrintStatus = VSFALSE ;
    Parm->compute_eig = VSTRUE ;
    
    /* seed for random number generator */
    Parm->seed = 1 ;

    /* machine epsilon for VSFLOAT */
    VSFLOAT eps, t ;
    eps = VSONE ;
    t = VSONE ;
    while ( t > 0 )
    {
        eps /= 2. ;
        t = VSONE + eps ;
        t -= VSONE ;
    }
    eps *= 2 ;
    Parm->eps = eps ;

    /* generate solution using full diagonalization of A when n <= SmallMatrix*/
    Parm->SmallMatrix = 16 ;

    /* When n > SmallMatrix, solution is generated using partial diagonalization
       based on the Lanczos process.  */
    Parm->basis_factor = 4 ;       /* # Lanczos iterations = basis_factor*
                                        log2(matrix dim) */
    Parm->Lanczos_tol = eps ;      /* |u_j| <= Lanczos_tol*Amax => stop Lanc. */
    Parm->qr_its = 4 ;             /* qr_its times dim = total # qr iterations*/
    Parm->qr_tol_rel = eps ;       /* relative eigenvalue convergence tol*/
    Parm->qr_tol_abs = eps*1.e-3  ;/* absolute eigenvalue convergence tol*/

    /* lower bound for qr scaling factor; if this is made large, then
       allocation for DT->gi and DT->gs should be increased so that
       -log_2(qr_lower) < the factor in these allocation */
    Parm->qr_lower = 1.e-9 ;

    /* number of gradient projection steps */
    Parm->ngp = 3 ;

    Parm->initial_project = VSFALSE ;
    Parm->try_both_covers = VSTRUE ;
}

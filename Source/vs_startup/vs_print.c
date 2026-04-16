/* ==========================================================================
   === vs_print_status ====================================================
   ==========================================================================
    Print the status at termination of the run
   ========================================================================== */

#include "VS.h"

void VS_print_status
(
    int status  /* status */
)
{
    printf ("\nVS run status (Version %d.%d.%d, %s):\n\n",
               VS_MAIN_VERSION, VS_SUB_VERSION, VS_SUBSUB_VERSION,
               VS_DATE) ;

    if ( status == VS_QR_DIAGONALIZATION_FAILS )
    {
        printf ("In vs, unable to diagonalize the matrix\n") ;
    }
    else if ( status == VS_NOT_ENOUGH_SPACE_FOR_QR_DIAGONALIZATION )
    {
        printf ("In vs, unable to diagonalize the matrix\n") ;
    }
    else if ( status == VS_OUT_OF_MEMORY )
    {
        printf ("Out of memory in vs\n") ;
    }
}

/* ==========================================================================
   === vs_print_parms =======================================================
   ==========================================================================
   print parameter values
   ========================================================================== */
void vs_print_parms
(
    VSparm *Parm /* parameter structure */
)
{
    /* numerical parameters */
    printf ("print level (0 = none, 3 = maximum) .......... PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("machine epsilon ..................................... eps: %e\n",
             Parm->eps) ;
    printf ("random number seed ................................. seed: %i\n",
             Parm->seed) ;
    printf ("SmallMatrix max dimension ...............................: %i\n",
             Parm->SmallMatrix) ;
    printf ("Stop Lanczos when |u| <= Lanczos_tol*Amax .... Lanczo_tol: %e\n",
             Parm->Lanczos_tol) ;
    printf ("# Lanczos iter= basis_factor*log2(dimension) basis_factor: %e\n",
             Parm->basis_factor) ;
    printf ("stop QR when u_{k-1}^2 <= qr_tol_rel*|d_k| ... qr_tol_rel: %e\n",
             Parm->qr_tol_rel) ;
    printf ("stop QR when |u_{k-1}| <= qr_abs_tol ......... qr_tol_abs: %e\n",
             Parm->qr_tol_abs) ;
    printf ("QR iterations <= qr_its * dimension .............. qr_its: %e\n",
             Parm->qr_its) ;
    printf ("lower bound for QR scaling factor .............. qr_lower: %e\n",
             Parm->qr_lower) ;
    printf ("PrintParms ..............................................: ") ;
             vs_print_TF (Parm->PrintParms) ;
    printf ("PrintStatus .............................................: ") ;
             vs_print_TF (Parm->PrintStatus) ;
}

void vs_print_TF 
(
    int TF /* TRUE or FALSE (or EMPTY if not True or False) */
)   
{   
    if ( TF == VSTRUE )
    {
        printf ("TRUE\n") ;
    }
    else if ( TF == VSFALSE )
    {
        printf ("FALSE\n") ;
    }
    else
    {
        printf ("EMPTY\n") ;
    }
}

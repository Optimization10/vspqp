#include "shared_api.h"
#include "limits.h"

void VSP_DumpNeighs(VSPINT i, VSPINT *Ai, VSPLONG *Ap)
{
    FILE *file = fopen("output/dump.txt", "a") ;
    fprintf(file, "\nNeighbors of %lld:\n", (long long) i) ;
    for(VSPLONG p = Ap[i] ; p < Ap[i+1] ; p++)
    {
        VSPINT neigh = Ai[p] ;
        if(p != Ap[i+1] -1)
        {
            fprintf(file, "%lld,  ", (long long) neigh) ;
        }
        else
        {
            fprintf(file, "%lld", (long long) neigh) ;
        }

    }
    fclose(file) ;
    return ;
}

void VSP_ResetDumpFile()
{
    FILE *file = fopen("output/dump.txt", "w") ;
    fclose(file) ;
    return ;
}

void VSP_check_matrix_is_symmetric(
    VSPINT           *Ai,
    VSPINT           *Ap,
    VSPINT           *Ae,
    VSPLONG nz,
    VSPINT             n,
    VSP_options* options
) ;

/* Check each state individually, as well as the cross term */
void VSP_check_states(VSP_state *X, VSP_state *Y, VSP_problem *problem)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSP_check_state(X, problem) ;
    VSP_check_state(Y, problem) ;

    /* Check the cross term */
    VSPFLOAT *x_f = VSP_point_to_array(X->point) ;
    VSPFLOAT *y_f = VSP_point_to_array(Y->point) ;
    VSPFLOAT *Ax_f = VSP_Az_struct_to_array(X->Az) ;
    VSPFLOAT *Ay_f = VSP_Az_struct_to_array(Y->Az) ;

    VSPFLOAT xAy_1 = 0 ;
    VSPFLOAT xAy_2 = 0 ;

    for(VSPINT i = 0 ; i < problem->Graph->n ; i++)
    {
        xAy_1 += x_f[i] * Ay_f[i] ;
        xAy_2 += y_f[i] * Ax_f[i] ;
    }

    if( VSP_equal(xAy_1, xAy_2) == false )
    {
        VSP_error("Cross terms don't agree", DBG_INFO) ;
    }

    /* Free memory */
    free(x_f) ;
    free(y_f) ;
    free(Ax_f) ;
    free(Ay_f) ;
}

void Check_Part
(
    VSPINT *x,
    VSPINT *Ax,
    VSPINT *y,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,
    VSP_problem *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;
    VSPINT n = problem->Graph->n ;
    /*------------------------------------------------------------------------*/
    // Load the data into states and use VSP_check_states
    /*------------------------------------------------------------------------*/
    VSP_point xPoint, yPoint ;

    xPoint.int_part = x ;
    xPoint.ind_frac = EMPTY ;
    xPoint.val_frac = 0 ;
    xPoint.n = n ;
    
    yPoint.int_part = y ;
    yPoint.ind_frac = EMPTY ;
    yPoint.val_frac = 0 ;
    yPoint.n = n ;
    
    VSP_Az Axstruct, Aystruct ;

    Axstruct.int_part    = Ax ;
    Axstruct.ind_frac_z  = EMPTY ;
    Axstruct.val_frac_z  = 0 ;
    Axstruct.n           = n ;
    Axstruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Axstruct.one         = 1 ;


    Aystruct.int_part    = Ay ;
    Aystruct.ind_frac_z  = EMPTY ;
    Aystruct.val_frac_z  = 0 ;
    Aystruct.n           = n ;
    Aystruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Aystruct.one         = 1 ;

    VSP_state X, Y ;

    X.point  = &xPoint ;
    X.Az     = &Axstruct ;
    X.weight = wgtx ;

    Y.point  = &yPoint ;
    Y.Az     = &Aystruct ;
    Y.weight = wgty ;

    VSP_check_states(&X, &Y, problem) ;


    if(wgtx < 0 || wgty < 0)
    {
        ERR("Negative weights") ;
    }

    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;
    return ;
}

void VSP_check_max_of_array
(
    VSPINT        *c,
    VSPINT         n, /* Length of c */
    VSPINT max_given /* Purported maximum value of c */
) ;


/* Checks:
    - The point is feasible
    - Az is correct
    - The weight is correct   */
void VSP_check_state(VSP_state *Z, VSP_problem *problem)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    /* Check that the n's agree */
    VSPINT nFromZ       = Z->point->n ;
    VSPINT nFromAz      = Z->Az->n ;
    VSPINT nFromProblem = problem->Graph->n ;
    if(nFromZ != nFromProblem)
    {
        VSP_error("n is incorrect", DBG_INFO) ;
    }

    if(nFromAz != nFromProblem)
    {
        VSP_error("n is incorrect", DBG_INFO) ;
    }

    /* Check the point */
    VSP_check_point(Z->point) ;

    /* Check that Az is correct */
    VSP_check_Az(Z->point, Z->Az, problem) ;

    /* Check that the weight is correct, and z satisfies constraint */
    VSP_check_weight(Z, problem->Graph->c) ;
    if(Z->weight > problem->u)
    {
        VSP_error("Upper balance violated ", DBG_INFO);
    }

    return ;
}



void VSP_check_weight(VSP_state *Z, VSPINT *c) {

    static long long unsigned count = 0 ;
    count++ ;

    /* Local vars */
    VSPINT i ;
    VSPINT n = Z->point->n ;
    VSPINT *z_int = Z->point->int_part ;
    VSPINT weight_check_int = 0 ;
    double   weight_check_frac = 0 ;

    /* Compute the weight of the integer part */
    for (i = 0 ; i < n ; i++)
    {
        weight_check_int += c[i] * z_int[i] ;
    }

    /* Compute the weight of the fractional part */
    VSPINT ind_frac = Z->point->ind_frac ;
    if(ind_frac != EMPTY) {
        weight_check_frac = c[ind_frac] * Z->point->val_frac ;
    }

    if( !VSP_equal(weight_check_int + weight_check_frac, Z->weight) )
    {
        VSP_error("Incorrect weight", DBG_INFO) ;
    }
}

void VSP_check_Az
(
    VSP_point         *z,
    VSP_Az           *Az,
    VSP_problem *problem
)
{

    static long long unsigned count = 0 ;
    count ++ ;

    VSPINT i ;
    VSPINT n       = problem->Graph->n ;
    VSPINT *Az_int = Az->int_part ;

    /* Compute Az from scratch */
    VSP_Az *Az_dbg = VSP_construct_Az(n) ;
    Compute_Az(z, Az_dbg, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;

    /* Check that each entry in the integer part is in [0, n] */
    for(i = 0 ; i < n ; i++) {
        if(Az_int[i] < 0 || Az_int[i] > n) {
            VSP_error("Problem in Az", DBG_INFO) ;
            break ;
        }
    }

    /* Check the integer parts are equal */
    VSP_check_arrays_are_equal(Az->int_part, Az_dbg->int_part, n) ;

    /* Check that the fractional part is equal */
    if(Az->ind_frac_z != Az_dbg->ind_frac_z) {
        VSP_error("Error in Az", DBG_INFO) ;
    }

    if(Az->ind_frac_z != EMPTY) {

        /* Check that the integer part of the fractional entry is < n*/
        if(Az_int[Az->ind_frac_z] > n - 1 || Az_int[Az->ind_frac_z] < 0) {
            VSP_error("Problem with fractional index in Az", DBG_INFO) ;
        }

        /* Check that the fractional columns are the same */
        bool    same              = true ;
        bool    one_is_ok         = true ;
        VSPINT *col_frac_z        = Az->col_frac_z ;
        VSPINT *col_frac_z_dbg    = Az_dbg->col_frac_z ;
        VSPINT  one               = Az->one ;
        VSPINT  one_dbg           = Az_dbg->one ;

        for(i = 0 ; i < n ; i++) {
            bool Az_is_one     = (col_frac_z[i]     == one) ;
            bool Az_dbg_is_one = (col_frac_z_dbg[i] == one_dbg) ;

            if(Az_is_one != Az_dbg_is_one) {
                same = false ;
            }
            /* Check that each entry is between 0 and one */
            if(col_frac_z[i] < 0 || col_frac_z[i] > one) one_is_ok = false ;
            if(col_frac_z_dbg[i] < 0 || col_frac_z_dbg[i] > one_dbg) {
                one_is_ok = false ;
            }
        }

        if(same == false)
        {
            VSP_error("Fractional columns mismatch", DBG_INFO) ;
        }
        if(one_is_ok == false) {
            VSP_error("Problem with fractional column", DBG_INFO) ;
        }

        /* Check that the fractional column has a one in index ind_frac_z */
        if(col_frac_z[Az->ind_frac_z] != one) {
            VSP_error("Matrix column wrong", DBG_INFO) ;
        }
    }

    /* Free memory */
    VSP_destruct_Az(Az_dbg) ;

    return ;
}



void VSP_check_arrays_are_equal
(
    const VSPINT *const u, /* Check that vector u equals vector v */
    const VSPINT *const v,
    const VSPINT  n  /* Size of u and v */
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i ;
    bool equal = true ;
    for(i = 0 ; i < n ; i++)
    {
        if(u[i] != v[i])
        {
            equal = false ;
            break ;
        }
    }
    if(equal == false)
    {
        VSP_error("Arrays not equal", DBG_INFO) ;
    }
    return ;
}

/*
 INPUTS:
 a) A matrix A in compressed column format
 b) Edge weights (optional)

 OUTPUTS:
 a) If the matrix passes all the checks, return TRUE (i.e. integer 1).
    Otherwise, return FALSE (integer 0) and print errors to the error log.
*/
VSPINT VSP_check_matrix
(
    VSPINT* Ai,
    VSPINT* Ap,
    VSPINT* Ae,
    VSPLONG nz,
    VSPINT n,
    VSP_options* options
) {
    static long long unsigned count = 0 ;
    count += 1;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, j, p, k, end ;
    VSPINT repeated_entry             = FALSE ;
    VSPINT diag_nonzero               = FALSE ;
    VSPINT invalid_neighbor           = FALSE ;
    VSPINT error                      = FALSE ;
    VSPINT invalid_entry_in_Ap        = FALSE ;
    VSPINT last_entry_of_Ap_incorrect = FALSE ;



    /* Check that Ap is (weakly) increasing, and each entry is in [0, nz]*/
    VSPINT Ap_not_increasing = FALSE ;
    for (i = 0 ; i < n ; i++) {
        if ( Ap [i] > Ap [i+1] ) Ap_not_increasing = TRUE ;
        if ( Ap [i] < 0 || Ap [i] > nz) invalid_entry_in_Ap = TRUE ;
    }

    /* Check that the last entry of Ap is nz */
    if ( Ap [n] != nz ) last_entry_of_Ap_incorrect = TRUE ;


    /* Check the following: Members of Ai are in [0, n-1], no repeated entries
    any chunk representing a column, diagonal is 0 */
    VSPINT* col = malloc(n * sizeof(VSPINT)) ;
    for (k = 0 ; k < n ; k++) {

        /* zero out the column*/
        for (j = 0 ; j < n ; j++) col [j] = 0 ;

        end = Ap [ k + 1 ] ;
        for (p = Ap [k] ; p < end ; p++) {
            i = Ai [p] ;

            /* Check that neighbors are not repeated */
            if( col [i] == 0 ) {
                col [i] = 1 ;
            }
            else {
                repeated_entry = TRUE ;
            }

            /* Check that the neighbor is between 0 and n-1, inclusive */
            if ( i < 0 || i > n - 1) {
                invalid_neighbor = TRUE ;
            }


            /* Check that the edge cost is at least 1, and reasonably small */
            if ( Ae [p] < 1 ) {
                VSP_error(" Error in edge weights", DBG_INFO) ;
            }
            if ( Ae [p] > 1E10 ) {
                VSP_error("Warning: probable error in edge weights", DBG_INFO) ;
            }

        }

        /* Check that the diagonal is 0 */
        if (col [k] != 0) {
            diag_nonzero = TRUE ;
        }
    }

    /* Check that the matrix is symmetric */
    VSP_check_matrix_is_symmetric(Ai, Ap, Ae, nz, n, options) ;

    /* Print errors, if there are any */
    if ( Ap_not_increasing ) {
        VSP_error("Ap not increasing", DBG_INFO) ;
    }

    if ( repeated_entry ) {
        VSP_error("Repeated neighbor in Ai", DBG_INFO) ;
        error = TRUE ;
    }

    if ( diag_nonzero ) {
        VSP_error("Nonzero diagonal", DBG_INFO) ;
        error = TRUE ;
    }

    if ( invalid_neighbor ) {
        VSP_error("Invalid neighbor in graph", DBG_INFO) ;
        error = TRUE ;
    }

    if ( last_entry_of_Ap_incorrect ) {
        VSP_error("Last entry of Ap incorrect", DBG_INFO) ;
    }

    if ( invalid_entry_in_Ap ) {
        VSP_error("Invalid entry in Ap", DBG_INFO) ;
    }

    /* Free memory */
    free(col) ;

    return error ;
}


void VSP_check_matrix_is_symmetric
(
    VSPINT           *Ai,
    VSPINT           *Ap,
    VSPINT           *Ae,
    VSPLONG nz,
    VSPINT             n,
    VSP_options* options
) {
    /* Get the transpose */
    VSP_problem* T = VSP_transpose(Ai, Ap, Ae, nz, n, options) ;
    VSPINT *Ti = T->Graph->Ai, *Te = T->Graph->Ae ;
    VSPLONG *Tp = T->Graph->Ap ;

    /* Check that the number of 1's in each column is the same in A and T */
    VSP_check_arrays_are_equal(Ap, Tp, n+1) ;

    /* Check that the row index of each 1 in each column is the same in A, T*/
    for(VSPINT j = 0 ; j < n ; j++) {
        VSP_check_permutation(&Ai[Ap[j]], &Ti[Tp[j]], Ap[j+1] - Ap[j]) ;
    }

    if (Ae != NULL) VSP_check_permutation(Ae, Te, nz) ;

    /* Free memory */
    VSP_destruct_problem(T) ;
}


/* ============================================================================
   ================ VSP_check_entries_are_between_bounds ======================
   ============================================================================

 PURPOSE:
 Given an array of ints and lower and upper bounds, check that every entry is
 between the bounds


 INPUTS:
 a) An array of ints
 b) Lower and upper bounds

 OUTPUTS:
 a) If there exists an entry which is not between the bounds, an error is
    written to the error log
*/

void VSP_check_entries_are_between_bounds(
    VSPINT*  v,
    VSPINT   n, /* Length of v */
    VSPINT  lb,
    VSPINT  up
) {
    WARN_IF_NOT_DEBUG_MODE ;
    static unsigned long long count = 0 ;
    count += 1 ;
    if (count == ULLONG_MAX) {
        printf("\n\n WARNING: Max iteration reached in %s \n\n", __func__) ;
        exit(1) ;
    }

    VSPINT i, vi ;
    VSPINT error = 0 ;

    /* Compute the true maximum of c */
    for (i = 0 ; i < n ; i++) {
        vi = v [i] ;
        if ( vi < lb || vi > up) {
            error = 1 ;
        }
    }

    if ( error ) {
        VSP_error("Array entry out of bounds", DBG_INFO) ;
    }

}

/* ============================================================================
   ================ VSP_check_max_of_array ====================================
   ============================================================================

 PURPOSE:
 Given an array of ints and its purported maximum entry, check that the given
 maximum is correct


 INPUTS:
 a) An array of ints
 b) A purported maximum entry

 OUTPUTS:
 a) If the given maximum is not correct, an error is written to the error log
*/

void VSP_check_max_of_array(
    VSPINT        *c,
    VSPINT         n, /* Length of c */
    VSPINT max_given /* Purported maximum value of c */
) {

    static unsigned long long count = 0 ;
    count += 1 ;

    WARN_IF_NOT_DEBUG_MODE ;

    if (count == ULLONG_MAX) {
        printf("\n\n WARNING: Max iteration reached in ") ;
        printf(__func__) ;
        printf("\n\n") ;
        exit(1) ;
    }

    VSPINT i ;
    VSPINT max_true = c [0] ;

    /* Compute the true maximum of c */
    for (i = 0 ; i < n ; i++) {
        if (c [i] > max_true) max_true = c [i] ;
    }

    if (max_given != max_true) {
        VSP_error("Max of array not correct", DBG_INFO) ;
    }
}

/* ============================================================================
   ================ VSP_check_problem =========================================
   ============================================================================

 PURPOSE:
 Check that all the data given corresponds to a valid VSP. In particular,
 check:
 - The given matrix is the adjacency matrix of an undirected Graph-> The edge
   weights, if they exist, are at least 1, and symmetric
 - The costs and weight arrays have entries >= 1
 - gamma is the maximum of the cost array
 - The x and y upper balance constraints are the same, and are at least 1
 - The sum of all the costs matches problem->totalgraphweight
 - The max degree stored in the problem is correct

 INPUTS:
 a) A VSP_problem struct

 OUTPUTS:
 a) If the matrix passes all the checks, return TRUE (i.e. integer 1).
    Otherwise, return FALSE (integer 0) and print errors to the error log.
*/

VSPINT VSP_check_problem
(
    VSP_problem* problem,
    VSP_options* options
)
{
    static long long unsigned count = 0 ;
    count += 1 ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, deg ;

    /* Unpack problem */
    VSPINT* Ai = problem->Graph->Ai ;
    VSPINT* Ap = problem->Graph->Ap ;
    VSPINT* Ae = problem->Graph->Ae ;
    VSPINT* c = problem->Graph->c ;

    VSPINT alpha = problem->Graph->min_vtx_cost ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;
    VSPLONG nz    = problem->Graph->nz ;
    VSPINT n     = problem->Graph->n ;
    VSPINT bal   = problem->u ;

    VSPINT cost_is_wrong     = FALSE ;
    VSPINT weight_is_wrong   = FALSE ;
    VSPINT is_there_an_error = FALSE ;


    /* Do all the matrix checks */
    VSP_check_matrix(Ai, Ap, Ae, nz, n, options) ;

    /* Check that ci and wi are >= 1. Also, compute alpha, gamma, and the
       sum of all ci, and max degree */
    VSPINT locally_computed_alpha = c [0] ;
    VSPINT locally_computed_gamma = EMPTY ;
    VSPINT sum_of_all_ci = 0 ;
    VSPINT maxDeg = EMPTY ;
    for (i = 0 ; i < n ; i++)
    {
        VSPINT ci = c [i] ;
        sum_of_all_ci += ci ;
        if (ci > locally_computed_gamma) locally_computed_gamma = ci ;
        if (ci < locally_computed_alpha) locally_computed_alpha = ci ;
        if (ci < 1) cost_is_wrong = TRUE ;

        deg = Ap[i+1] - Ap[i] ;
        if(deg > maxDeg)
        {
            maxDeg = deg ;
        }
    }

    if(maxDeg != problem->maxDegree)
    {
        ERR("Max degree of problem is wrong") ;
    }

    if (cost_is_wrong) {
        VSP_error("Cost is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (weight_is_wrong) {
        VSP_error("Weight is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (locally_computed_alpha != alpha) {
        VSP_error("Gamma is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (locally_computed_gamma != gamma) {
        VSP_error("Gamma is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (bal < 1) {
        VSP_error("Balance constraints are wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (sum_of_all_ci != problem->Graph->total_vtx_cost) {
        VSP_error("Problem with sum of all ci", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }
    if(problem->ratioMaxCostNodes <= 0 || problem->ratioMaxCostNodes > 1)
    {
        VSP_error("Wrong ratio of max weight nodes", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    return is_there_an_error ;

}



/* Check whether arr_2 is a permutation of arr_1 */
void VSP_check_permutation (
    VSPINT* arr_1,
    VSPINT* arr_2,
    VSPINT n
) {
    static long long unsigned count = 0 ;
    count += 1 ;


    VSPINT i ;

    /* Make copies of the input arrays */
    VSPINT* arr_1_copy = malloc (n * sizeof(VSPINT)) ;
    VSPINT* arr_2_copy = malloc (n * sizeof(VSPINT)) ;

    for (i = 0 ; i < n ; i++) {
        arr_1_copy [i] = arr_1 [i] ;
        arr_2_copy [i] = arr_2 [i] ;
    }

    /* Sort the copies */
    qsort(arr_1_copy, n, sizeof(*arr_1), VSP_comparison_increasing) ;
    qsort(arr_2_copy, n, sizeof(*arr_2), VSP_comparison_increasing) ;

    VSPINT arrays_are_equal = TRUE ;
    for (i = 0 ; i < n ; i++) {
        if (arr_1_copy [i] != arr_2_copy [i]) {
            arrays_are_equal = FALSE ;
        }
    }

    if (arrays_are_equal == FALSE) {
        VSP_error("Array not a permutation", DBG_INFO ) ;
    }

    /* Free memory */
    free(arr_1_copy) ;
    free(arr_2_copy) ;
    return ;
}


/* Check as much as possible when only the Az struct is available             */
void VSP_light_check_Az(
    VSP_Az *Az
) {

    static long long unsigned count = 0 ;
    count++ ;

    VSPINT  n      = Az->n ;
    VSPINT *Az_int = Az->int_part ;

    /* Check that the integer part is between 0 and n */
    for(VSPINT i = 0 ; i < n ; i++) {
        if(Az_int[i] < 0 || Az_int[i] > n) {
            VSP_error("Error in Az", DBG_INFO) ;
        }
    }


    /* Check that the fractional column has a one in index ind_frac_z */
    if(Az->ind_frac_z != EMPTY)
        {
        if(Az->col_frac_z[Az->ind_frac_z] != Az->one) {
            VSP_error("Matrix column wrong", DBG_INFO) ;
        }
    }
}


void VSP_check_point(VSP_point *z) {

    static long long unsigned count = 0 ;
    count++ ;

    VSPINT *z_int = z->int_part ;
    VSPINT  n     = z->n ;

    /* Check that the integer part is binary */
    bool binary = true ;
    for(VSPINT i = 0 ; i < n ; i++) {
        if (z_int[i] != 0 && z_int[i] != 1) {
            binary = false ;
        }
    }
    if( binary == false ) VSP_error("Not binary", DBG_INFO) ;


    /* Check that the integer part of the fractional entry is 0 */
    VSPINT ind_frac = z->ind_frac ;
    if(ind_frac != EMPTY && z_int[ind_frac] != 0) {
        VSP_error("Problem with integer part", DBG_INFO) ;
    }

    /* Check that the fractional index is empty or in the valid range */

    if( ind_frac != EMPTY && !(0 <= ind_frac && ind_frac <= n-1) ) {
        VSP_error("Fractional index out of bounds", DBG_INFO) ;
    }

    /* Check that the fractional value is in [0, 1] */
    VSPFLOAT val_frac = z->val_frac ;
    if(val_frac < 0 || val_frac > 1) {
        VSP_error("Issue with fractional value", DBG_INFO) ;
    }

    /* If frac index exists, check that the frac value is in (0, 1)*/
    if(ind_frac != EMPTY) {
        if(val_frac <= 0 || val_frac >= 1 ) {
            VSP_error("Issue with fractional value", DBG_INFO) ;
        }
    }
    return ;
}

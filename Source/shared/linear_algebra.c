#include "core.h"
#include "shared_api.h"
#include "shared_functions.h"
#include "shared_macros.h"

void VSP_GetSepCost
(
    VSP_state          *X,
    VSP_state          *Y,
    Problem  *problem,
    VSPINT       *sepCost   /* Output*/
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT  n = X->point->n ;
    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;
    VSPINT *c = problem->Graph->c ;

    VSPINT localSepCost = 0 ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        if(x[i] == 0 && y[i] == 0)
        {
            localSepCost += c[i] ;
        }
    }
    *sepCost = localSepCost ;
    return ;
}

/* Given a binary point, return xAy */
void VSP_num_pen
(
    VSP_state *X,
    VSP_state *Y,
    Problem *problem,
    VSPFLOAT *numPen
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT *x = X->point->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT xAy = 0 ;

    


    for(i = 0 ; i < n ; i++)
    {
        xAy += x[i] * Ay[i] ;
    }

    #ifdef VSP_DEBUG
        VSPINT yAx = 0 ;
        for(i = 0 ; i < n ; i++)
        {
            yAx += Y->point->int_part[i] * X->Az->int_part[i] ;
        }
        if(xAy != yAx)
        {
            Error("Error in mixed term", DBG_INFO) ;
        }
    #endif

    *numPen = xAy ;
    return ;
}


/* Given a binary point, return xAy */
void VSP_xAy
(
    VSP_state         *X,
    VSP_state         *Y,
    Problem *problem,
    VSPFLOAT        *xAy
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPFLOAT localxAy, xiFrac, AyiFrac, xi, Ayi ;

    VSPINT *xInt = X->point->int_part ;
    VSPINT *AyiInt = Y->Az->int_part ;
    
    VSPINT xHasFrac = (X->point->ind_frac != EMPTY) ;
    VSPINT AyHasFrac = (Y->Az->ind_frac_z != EMPTY) ;

    localxAy = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        if(xHasFrac)
        {
            xiFrac = (X->point->ind_frac == i) ? X->point->val_frac : 0 ;
        }
        else
        {
            xiFrac = 0 ;
        }

        xi = xInt[i] + xiFrac ;

        if(AyHasFrac)
        {
            if( Y->Az->col_frac_z[i] == Y->Az->one )
            {
                AyiFrac = Y->point->val_frac ;
            }
            else
            {
                AyiFrac = 0 ;
            }
        }
        else
        {
            AyiFrac = 0 ;
        }

        Ayi = AyiInt[i] + AyiFrac ;
        localxAy += xi * Ayi ;
    }
    *xAy = localxAy ;
    return ;
}

/* Compute c - gamma Aw and store it in grad_z. No assumption on contents
   of grad_z */
void VSP_compute_gradient_v1
(
    VSPFLOAT *grad_z,
    VSP_Az       *Aw,
    VSPINT        *c,
    VSPFLOAT   gamma
)
{
    VSPINT i ;
    VSPINT n = Aw->n ;
    VSPINT *Aw_int = Aw->int_part ;

    if(Aw->ind_frac_z != EMPTY)
    {
        VSPINT   one      = Aw->one ;
        VSPINT  *col_frac = Aw->col_frac_z ;
        VSPFLOAT val_frac = Aw->val_frac_z ;
        for(i = 0 ; i < n ; i++)
        {
            VSPFLOAT Awi = Aw_int[i] + val_frac * ( col_frac[i] == one ) ;
            grad_z[i] = c[i] - gamma * Awi ;
        
        }
    }

    for(i = 0 ; i < n ; i++)
    {
        grad_z[i] = c[i] - gamma * Aw_int[i] ;
    }

    return ;
}



/* On input:
   Compute (A + I)*z where z is a fully binary vector
   1) No assumptions on contents of z or Az arrays

   On output:
   2) The vector Az is (A + I) * z */


// #define USE_MKL
void VSP_mat_vec_mult
(
    const VSPINT *const    z,
    VSPINT       *const   Az,
    const VSPINT *const   Ai,
    const VSPLONG *const  Ap,
    const VSPINT           n
)
{
    TIMER_START ;


    VSPINT i, col = 0 ;
    VSPLONG p ;
    VSPLONG start, end ;

    /* Initialize Az to all 0's */
    for(i = 0; i < n; i++) Az[i] = 0 ;

    for (col = 0; col < n; col++)
    {
        if (z[col])
        {
            start = Ap[col] ;
            end   = Ap[col + 1] ;
            for (p = start ; p < end ; p++) Az[ Ai[p] ] ++ ;
            /* Take identity into account */
            Az[col] ++ ;
        }
    }
    TIMER_END ;
}

/* Update Az struct, given z as a point struct */
void Compute_Az
(
    VSP_point              *z,
    VSP_Az                *Az,
    VSPINT                *Ai,
    VSPLONG               *Ap,
    VSPINT                  n
)
{
    SET_FUNCTION_COUNTER ;

    // Compute integer part
    VSP_mat_vec_mult(z->int_part, Az->int_part, Ai, Ap, n) ;
    Az->ind_frac_z = z->ind_frac ;

    // Compute fractional part
    Az->ind_frac_z = z->ind_frac ;
    Az->val_frac_z = z->val_frac ;

    if (z->ind_frac != EMPTY)
    {
        Az->one++ ;
        VSP_extract_matrix_col(Az->col_frac_z, z->ind_frac, Az->one, Ai, Ap) ;
    }
}

// void Az_Compute
// (
//     VSPINT     *zint,
//     VSPINT   indFrac,
//     VSPFLOAT valFrac,
//     VSP_Az       *Az,
//     VSPINT       *Ai,
//     VSPLONG      *Ap,
//     VSPINT         n
// )
// {
//     SET_FUNCTION_COUNTER ;
    
//     // Compute the integer part of Az
//     VSP_mat_vec_mult(zint, Az->int_part, Ai, Ap, n) ;

//     if(indFrac != EMPTY)
//     {
//         Az->ind_frac_z = indFrac ;
//         Az->val_frac_z = valFrac ;
//         Az->one++ ;
//         // VSP_extract_matrix_col(Az->col_frac_z, z->ind_frac, Az->one, Ai, Ap) ;
//     }
// }


void VSP_Az_pm_col_of_A (
    VSPINT             *Az, // vector to which we add/subtract
    VSPINT               i, // column to add/subtract
    VSPINT          operation, // 1 if add, -1 if subtract
    Problem   *problem
) {
    VSPINT *Ai = problem->Graph->Ai ;
    VSPLONG *Ap = problem->Graph->Ap ;
    VSPLONG const end = Ap [i+1] ;
    for(VSPLONG p = Ap[i] ; p < end ; p++)
    {
        VSPINT const j = Ai [p] ;
        Az[j] += operation;
    }
    // don't forget the identity matrix
    Az[i] += operation;
}



/* Compute the weight of a point. This function assumes that the true weight
   is an integer.                                                             */
void VSP_compute_weight(VSP_point *Z, VSPINT *weight_z, VSPINT *c) {
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    VSPINT  i ;
    VSPINT  weight = 0 ;
    VSPINT *z_int  = Z->int_part ;
    VSPINT  n      = Z->n ;

    /* Compute the weight of the integer part */
    for(i = 0 ; i < n ; i++) {
        weight += z_int[i] * c[i] ;
    }

    /* Add the weight of the fractional part */
    if(Z->ind_frac != EMPTY) {
        VSPFLOAT weight_frac_part = Z->val_frac * c[Z->ind_frac] ;
        /* TODO safety: Consider safety here. Maybe introduce a sort of check
                 before typecasting */
        weight += lroundf(weight_frac_part) ;
        #ifdef VSP_DEBUG
            if( fabs(weight_frac_part - round(weight_frac_part)) > VSP_ZERO ) {
                Error("Error with weight", DBG_INFO) ;
            }
        #endif
    }
    *weight_z = weight ;
}

Problem* VSP_transpose(
    VSPINT           *Ai, /* Row indices of non-zeros  */
    VSPLONG *Ap, /* Column pointers           */
    VSPINT           *Ae, /* Entries                   */
    VSPLONG nz, /* Number of rows/columns    */
    VSPINT             n, /* Number of nonzeros        */
    VSP_options* options
) {
 
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    #ifdef VSP_DEBUG
        if (Ap [n] != nz) {
            Error("Nz doesn't match Ap[n]", __FILE__, __LINE__, count) ;
        }
    #endif

    VSPLONG i ;
    VSPINT j, ind_row ;
    VSPLONG p ;
    VSPLONG* pointers_into_Ti ; 
    VSPLONG* num_nonzeros_in_row = calloc(n, sizeof(*num_nonzeros_in_row)) ;

    /* Prepare the transpose matrix T */
    Problem* T = VSP_construct_problem(n, nz, NULL, NULL, NULL, options) ;
    VSPINT *Ti = T->Graph->Ai , *Te = T->Graph->Ae ;
    VSPLONG *Tp = T->Graph->Ap ;

    /* Count the number of nonzeros in each row */
    for (VSPLONG ii = 0 ; ii < nz ; ii++) {
        ind_row = Ai [ii] ;
        num_nonzeros_in_row [ ind_row ] += 1 ;
    }


    /* Set Tp[0] = 0 and Tp[i] = number of nonzeros in rows 0, 1, ..., i-1 */
    /* Then, set num_nonzeros_in_row [i] = Tp [i] for i = 0, ..., n-1 */
    VSP_partial_sums_of_array(num_nonzeros_in_row, Tp, n) ;

    /* Change the name of the array, since the array has changed.
 
       Recall that Ti is split into chunks, corresponding to columns.
       pointers_in_Ti [i] holds the index of the first available entry
       in chunk i. Thus, initially, pointers_into_Ti [i] = Tp [i] for i = 
       0, 1, ..., n-1, and at the end, pointers_into_Ti [i] should
       equal Tp [i+1] */
    pointers_into_Ti = num_nonzeros_in_row ;

    /* Go through the entire matrix A, filling out the Ti and Te arrays */
    for ( j = 0 ; j < n ; j++) {
        for (p = Ap [j] ; p < Ap [j+1] ; p++) {
            i = Ai [p] ;
            
            /* A(i, j) = nonzero. Set T(j, i) = nonzero.
               Recall that Ae [p] is the entry in A(i, j). This value should
               become the entry in T(j, i) */
            Ti [ pointers_into_Ti [i] ] = j ;

            if (Ae != NULL ) {
                Te [ pointers_into_Ti [i] ] = Ae [p] ;
            }
            pointers_into_Ti [i] ++ ;
        }
    }

    #ifdef VSP_DEBUG

        VSPINT error_in_transpose = FALSE ;
        for (i = 0 ; i < n ; i++ ) {
            if ( pointers_into_Ti [i] != Tp [i + 1] ) {
                error_in_transpose = TRUE ;
            }
        }
        if (error_in_transpose) {
            Error("Error in transpose", DBG_INFO) ;
        }

        /* If the matrix is small enough, expand it out and check transpose */
        if (n < 300) {

            /* Create the explicit matrices */
            VSPINT** A_explicit = VSP_construct_explicit_matrix(NULL, Ai, nz,
                                                                      Ap, n  ) ;
            VSPINT** T_explicit = VSP_construct_explicit_matrix(NULL, Ti, nz,
                                                                      Tp, n  ) ;

            VSPINT** A_vals_explicit = NULL ;
            VSPINT** T_vals_explicit = NULL ;

            if (Ae != NULL) {
                A_vals_explicit = VSP_construct_explicit_matrix(
                                                            Ae, Ai, nz, Ap, n) ;

                T_vals_explicit = VSP_construct_explicit_matrix(
                                                            Te, Ti, nz, Tp, n) ;
            }

            /* Check the transpose matches */
            VSPINT transpose_wrong = FALSE ;
            for (i = 0 ; i < n ; i++) {
                for (j = 0 ; j < n ; j++) {

                    if ( A_explicit[i][j] != T_explicit[j][i] ) {
                        transpose_wrong = TRUE ;
                    }

                    if (Ae != NULL) {
                        if ( A_vals_explicit[i][j] != T_vals_explicit[j][i] ) {
                            transpose_wrong = TRUE ;
                    }
                    }

                }
            }
            if (transpose_wrong) {
                Error("Transpose wrong", DBG_INFO) ;
            }

            /* Free memory */
            VSP_destruct_explicit_matrix(A_explicit, n) ;
            VSP_destruct_explicit_matrix(T_explicit, n) ;

            if (Ae != NULL) {
                VSP_destruct_explicit_matrix(A_vals_explicit, n) ;
                VSP_destruct_explicit_matrix(T_vals_explicit, n) ;
            }

        }

    #endif

    /* Free memory */
    free(num_nonzeros_in_row) ;

    return T ;
}

/* WARNING: This function modifies both of the input arrays.
   On output, partial_sums[0] = 0, and for i = 1, 2, ..., n,

        partial_sums[i] = arr[0] + arr[1] + ... + arr[i-1].

   The array arr is also overwritten, so that for i = 0, 1, 2, ..., n-1, 
        arr[i] = partial_sums[i].
   */
VSPLONG VSP_partial_sums_of_array(
    VSPLONG          *arr,
    VSPLONG *partial_sums,
    VSPINT              n
) {
    VSPINT i ;
    VSPLONG sum = 0 ;

    for (i = 0 ; i < n ; i++) {
        partial_sums [i] = sum ;
        sum += arr [i] ;

        arr [i] = partial_sums [i] ;
    }
    partial_sums [n] = sum ;

    return sum ;
}



/* This function computes
        c' (X + Y) - gamma * X' A Y
   efficiently.

   The computation is done in the following way:
   Write X = x + fx e_p
         Y = y + fy e_q
   where x, y are the integer parts, p, q are the indices
   of the fractional entries, and fx, fy are the values of the 
   fractional entries in x and y. The objective becomes
      f(x, y) = c'(X+Y) - gamma( x' A y + fy x' A_q + fx y' A_p + fx fy A_pq )

   The quantities c'X and c'Y are the weights of the states, which are already
   known. Every other term is computed separately, and all terms are summed
   at the end.
*/
VSPFLOAT VSP_compute_obj(
    VSP_state         *X,
    VSP_state         *Y,
    Problem *problem,
    VSPFLOAT       gamma
) {
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        /* Check input */
        VSP_check_states(X, Y, problem) ;
    #endif

    /* Local vars */
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT xAy, yA_p, xA_q, A_pq ;
    VSPFLOAT obj_val = 0 ;

    VSPINT   *x      = X->point->int_part ;
    VSPINT    p      = X->point->ind_frac ;
    VSPFLOAT  fx     = X->point->val_frac ;
    VSPINT   *A_p    = X->Az->col_frac_z ;
    VSPINT    one_p  = X->Az->one ;

    VSPINT   *y      = Y->point->int_part ;
    VSPINT    q      = Y->point->ind_frac ;
    VSPFLOAT  fy     = Y->point->val_frac ;
    VSPINT   *A_q    = Y->Az->col_frac_z ;
    VSPINT    one_q  = Y->Az->one ;
    VSPINT   *Ay     = Y->Az->int_part ;

    /*TODO speed: Fastest way would be to check the existence of the fractional
                  indices at the very beginning and treat each case separately*/
    xAy  = 0 ;
    yA_p = 0 ;
    xA_q = 0 ;
    for(i = 0 ; i < n ; i++) {
        xAy += x[i] * Ay[i] ;
        if(p != EMPTY) yA_p += y[i] * (A_p[i] == one_p) ;
        if(q != EMPTY) xA_q += x[i] * (A_q[i] == one_q) ;
    }
    obj_val = X->weight + Y->weight - gamma * ( xAy + fy * xA_q + fx * yA_p ) ;

    /* Take into account the cross term */
    if(p != EMPTY && q != EMPTY) {
        A_pq = (A_p[q] == one_p) ;
        obj_val -= gamma * fx * fy * A_pq ;
    }

    #ifdef VSP_DEBUG
        if(n < 200) {
            VSPFLOAT obj_val_dbg ;
            obj_val_dbg = VSP_compute_objective_0(X, Y, problem, gamma) ;
            if( fabs(obj_val - obj_val_dbg) > VSP_ZERO )
            {
                Error("Error in computing objective", DBG_INFO) ;
            }
        }
    #endif

    return obj_val ;
}


/* This function computes
         c' (X + Y) - gamma * X' A Y
   by uncompressing the matrix A and converting X, Y to arrays
   of floats                                                     */
VSPFLOAT VSP_compute_objective_0(
    VSP_state *X,
    VSP_state *Y,
    Problem *problem,
    VSPFLOAT gamma
) {
    /* Local vars */
    VSPINT i, j ;
    VSPLONG start, end, p ;
    VSPFLOAT  obj_val = 0 ;
    VSPINT    n       = problem->Graph->n ;
    VSPINT   *Ai      = problem->Graph->Ai ;
    VSPLONG *Ap      = problem->Graph->Ap ;
    VSPINT   *c       = problem->Graph->c ;


    /* Convert X, Y to arrays of floats */
    VSPFLOAT *x_arr = VSP_point_to_array(X->point) ;
    VSPFLOAT *y_arr = VSP_point_to_array(Y->point) ;

    /* Allocate memory for A */
    VSPINT** A = malloc( n * sizeof(VSPINT*) ) ;
    for (i = 0 ; i < n ; i++) {
        A[i] = calloc( n, sizeof(VSPINT) ) ; // Use calloc to 0-initialize
    }

    /* Fill out A with 1's from Ai and Ap */
    for (i = 0 ; i < n ; i++) {
        start = Ap [i] ;
        end = Ap [ i + 1 ] ;
        for (p = start ; p < end ; p++) {
            j = Ai [p] ;
            A [i][j] = 1;
        }
    }
    /* Set the diagonal to 1's */
    for (i = 0 ; i < n ; i++) A [i][i] = 1 ;

    /* Compute Ay */
    VSPFLOAT *Ay = calloc(n, sizeof(*Ay)) ;

    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            Ay [i] += A [i][j] * y_arr [j] ;
        }
    }

    /* Compute objective */
    obj_val = 0 ;

    for (i = 0 ; i < n ; i++) {
        obj_val += c [i] * ( x_arr[i] + y_arr[i] ) - gamma * x_arr [i] * Ay [i] ;
    }

    /* Free memory */
    free(x_arr) ;
    free(y_arr) ;
    free(Ay) ;
    for (i = 0 ; i < n ; i++) free( A [i] ) ;
    free (A) ;

    return obj_val ;
}

// x and y are binary vectors
// Intended for use in debugging, as normally I would pass c^Tx and c^Ty
// instead of computing them
void Compute_ObjValDbg
(
    VSPINT               *x,
    VSPINT              *Ax,
    VSPINT               *y,
    VSPINT              *Ay,
    Problem    *problem,
    VSPINT       *retObjVal
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, objval ;
    VSPINT n     = problem->Graph->n ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;
    VSPINT *c = problem->Graph->c ;

    objval = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        objval += ( c[i] * (x[i] + y[i]) - gamma * x[i] * Ay[i]) ;
    }

    VSPINT objval2 = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        objval2 += ( c[i] * (x[i] + y[i]) - gamma * y[i] * Ax[i]) ;
    }
    if(objval != objval2)
    {
        ERR("Problem when computing objective") ;
    }

    
    *retObjVal = objval ;
    return ;
}

// x and y are binary vectors
void Compute_ObjValSimple
(
    VSPINT               *x,
    VSPINT              *Ax,
    VSPINT             wgtx,
    VSPINT               *y,
    VSPINT              *Ay,
    VSPINT             wgty,
    Problem    *problem,
    VSPINT       *retObjVal
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, objval ;
    VSPINT n     = problem->Graph->n ;
    VSPINT gamma = problem->Graph->max_vtx_cost ;
    VSPINT *c    = problem->Graph->c ;

    objval = 0 ;


    
    for(i = 0 ; i < n ; i++)
    {
        objval += ( c[i] * (x[i] + y[i]) - gamma * x[i] * Ay[i] ) ;
    }

    #ifdef VSP_DEBUG
        VSPINT objval2 = 0 ;
        for(i = 0 ; i < n ; i++)
        {
            objval2 += ( c[i] * (x[i] + y[i]) - gamma * y[i] * Ax[i]) ;
        }
        if(objval != objval2)
        {
            ERR("Problem when computing objective") ;
        }
    #endif
    
    *retObjVal = objval ;
    return ;
}

/* Compute (c - gamma * Aw)' z */
VSPFLOAT VSP_compute_LP_objective(
    VSP_point    *z,
    VSPINT       *c,
    VSPFLOAT  gamma,
    VSP_Az      *Aw
) {
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPFLOAT obj            = 0 ;  
    VSPFLOAT val_frac_w     = Aw->val_frac_z ;
    VSPFLOAT Awi ;
    VSPINT   n              = z->n ;
    VSPINT   w_has_ind_frac = (Aw->ind_frac_z != EMPTY) ;
    VSPINT  *col_A          = Aw->col_frac_z ;
    VSPINT  *Aw_int         = Aw->int_part ;
    VSPINT  *z_int          = z->int_part ;
    VSPINT   one            = Aw->one ;

    for(VSPINT i = 0 ; i < n ; i++) {
        // Awi = Aw_int[i] + w_has_ind_frac * (col_A[i] == one) * val_frac_w ;
        Awi = VSP_Azi(i, Aw) ;
        obj += z_int[i] * (c[i] - gamma * Awi) ;
    }

    /* Take into account the fractional index of z */
    if(z->ind_frac != EMPTY) {
        VSPINT k = z->ind_frac ;
        VSPFLOAT Awk = Aw_int[k] + w_has_ind_frac*(col_A[k] == one)*val_frac_w ;
        obj += z->val_frac * (c[k] - gamma * Awk) ;
    }

    return obj ;

}

/* Compute c - gamma Aw where w is an array of floats */
void VSP_ComputeGradFlt
(
    VSPFLOAT  *grad,
    VSPINT       *c,
    VSPFLOAT  gamma,
    VSPFLOAT     *z,
    VSPINT      *Ai,
    VSPLONG *Ap,
    VSPINT        n
)
{
    VSPINT i, col = 0 ;
    VSPLONG p ;
    VSPLONG start, end ;

    /* Initialize Az to all 0's */
    for(i = 0; i < n; i++) grad[i] = 0 ;

    for (col = 0; col < n; col++)
    {
        grad[col] += c[col] ;
        if (z[col])
        {
            start = Ap[col] ;
            end   = Ap[col + 1] ;
            for (p = start ; p < end ; p++) grad[ Ai[p] ] -= (gamma * z[col]) ;
            /* Take identity into account */
            grad[col] += (gamma * z[col]) ;
        }
    }
    return ;
}
#include "core.h"
#include "shared_api.h"
#include "ds_api.h"
#include <string.h>

/* ************* Explicit Matrix ******************************************* **/
VSPINT** VSP_construct_explicit_matrix
(
    VSPINT *values,
    VSPINT     *Ai,
    UNUSED VSPLONG nz,
    VSPLONG *Ap,
    VSPINT       n
)
{
    static long long unsigned count = 0 ;
    count += 1 ;

    VSPINT i, j ;
    VSPLONG end, p ;

    /* Make pointers for the rows of A */
    VSPINT** A = malloc( n * sizeof(VSPINT*) ) ;
    if ( A == NULL ) {
        VSP_error("Memory allocation failed", DBG_INFO) ;
        exit(1) ;
    }

    /* Allocate space for each row */
    for (i = 0 ; i < n ; i++) {
        A [i] = malloc( n * sizeof(VSPINT) ) ;
        if ( A [i] == NULL ) {
            VSP_error("Memory allocation failed", DBG_INFO) ;
            exit(1) ;
        }
    }

    /* Zero out the matrix */
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) A [i][j] = 0 ;
    }

    /* Fill out the matrix with 1's from Ai and Ap */
    for (j = 0 ; j < n ; j++) {
        end = Ap [ j + 1 ] ;
        for (p = Ap[j] ; p < end ; p++) {
            i = Ai [p] ;
            if (values == NULL) A [i][j] = 1;
            else A [i][j] = values [p] ;
        }
    }

    /* Check that the matrix is symmetric */
    VSPINT is_symmetric = TRUE ;
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            if ( A[i][j] != A [j][i] ) {
                is_symmetric = FALSE ;
                printf("\nA[%lld][%lld] = %lld   while   A[%lld][%lld] = %lld", (long long) i, (long long) j, (long long) A[i][j], (long long) j, (long long) i, (long long) A[j][i]) ;
                fflush(stdout) ;
            }
        }
    }

    if (is_symmetric == FALSE)
    {
        if (values == NULL)
        {
            VSP_error("Adjacency matrix not symmetric", DBG_INFO) ;
        }
        else
        {
            VSP_error("Edge weights matrix not symmetric", DBG_INFO) ;
        }
    }

    return A ;
}

void VSP_destruct_explicit_matrix( VSPINT** A, VSPINT n ) {
    VSPINT i ;

    for ( i = 0 ; i < n ; i++ ) free (A [i]) ;
    free(A) ;

    return ;
}
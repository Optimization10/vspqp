#include "core.h"
#include "shared_api.h"
#include "data_structures.h"
#include <stdint.h>
#include <string.h>
#include "utils.h"
#include <limits.h>


bool Is_Coarsest_Problem(Problem *problem)
{
    return (problem->child_problem == NULL) ; 
}

bool Is_Finest_Problem(Problem *problem)
{
    return (problem->parent_problem == NULL) ; 
}

VSPINT Are_Neighbors
(
    VSPINT vert1,
    VSPINT vert2,
    VSPINT   *Ai,
    VSPLONG *Ap
)
{
    for(VSPLONG p = Ap[vert1] ; p < Ap[vert1 + 1] ; p++)
    {
        if(Ai[p] == vert2) return 1 ;
    }
    return 0 ;
}

VSPINT penalty_term
(
    VSPINT  *z,
    VSPINT *Aw,
    VSPINT   n
)
{
    VSPINT ret = 0 ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        VSPINT term = z[i] * Aw[i] ;
        ret += term ;
    }
    return ret ;
}

VSPFLOAT zAw
(
    State *X,
    State *Y
)
{
    /* Local vars */
    VSPINT i, n = X->point->n, yA_p = 0, xA_q = 0, A_pq ;
    VSPFLOAT xAy = 0 ;
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
    if(p != EMPTY && q != EMPTY)
    {
        for(i = 0 ; i < n ; i++)
        {
            xAy  += x[i] * Ay[i] ;
            yA_p += y[i] * (A_p[i] == one_p) ;
            xA_q += x[i] * (A_q[i] == one_q) ;
        }
        xAy += fy * xA_q + fx * yA_p ;

        A_pq = (A_p[q] == one_p) ;
        xAy +=fx * fy * A_pq ;
    }
    else if(p != EMPTY)
    {
        for(i = 0 ; i < n ; i++)
        {
            xAy  += x[i] * Ay[i] ;
            yA_p += y[i] * (A_p[i] == one_p) ;
        }
        xAy += fx * yA_p ;
    }
    else if(q != EMPTY)
    {
        for(i = 0 ; i < n ; i++)
        {
            xAy += x[i] * Ay[i] ;
            xA_q += x[i] * (A_q[i] == one_q) ;
        }
        xAy += fy * xA_q ;
    }
    else
    {
        for(i = 0 ; i < n ; i++)
        {
            xAy += x[i] * Ay[i] ;
        }
    }

    return xAy ;
}

/* Count the number of nonzeros in a given array of integers */
VSPINT count_nzs(VSPINT *array, VSPINT n)
{
    VSPINT ret = 0 ;
    for(VSPINT i = 0 ; i < n ; i++)
    {
        if(array[i]) ret++ ;
    }
    return ret ;
}

void copy_array_int
(
    VSPINT *from,
    VSPINT *to,
    VSPINT  n
)
{
    VSPINT i ;
    for(i = 0 ; i < n ; i++) {
        to[i] = from[i] ;
    }
    return ;
}


/* Return the number of places in which a and b are different */
VSPINT bin_points_diff
(
    VSPINT *a,
    VSPINT *b,
    VSPINT n
)
{
    VSPINT i ;
    VSPINT n_diff = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        if(a[i] != b[i])
        {
            n_diff++ ;
        }
    }
    return n_diff ;

}

// Check if candidate is in array[start, ... end)
VSPINT Is_Array_Member
(
    VSPINT candidate,
    VSPINT    *array,
    VSPINT     start,
    VSPINT       end
)
{
    VSPINT i ;

    for (i = start ; i < end ; i++)
    {
        if (candidate == array [i]) return 1 ;
    }
    return 0 ;
}


bool Is_Complete_Graph
(
    VSPINT n, VSPLONG nz
)
{
    // Input validation
    if(n <= 0)
    {
        fprintf(stderr, "n <= 0 \n") ;
        return true ;
    }
    if(nz < 0)
    {
        fprintf(stderr, "nz < 0 \n") ;
        return true ;
    }

    long long nlong       = (long long) n ;
    long long nzlong      = (long long) nz ;
    long long upperbnd_nz = nlong * nlong - nlong ;

    if( nzlong ==  upperbnd_nz )
    {
        return true ;
    }
    else if( nzlong > upperbnd_nz )
    {
        fprintf(stderr, "n > n^2 - n\n") ;
        return true ;
    }
    else
    {
        return false ;
    }
}



/* Check if a == b */
bool equal
(
    VSPFLOAT a,
    VSPFLOAT b
)
{
    VSPFLOAT diff = fabs(a - b) ;
    if( diff < 1E-9 )
    {
        return true ;
    }

    VSPFLOAT A = fabs(a) ;
    VSPFLOAT B = fabs(b) ;
    VSPFLOAT largest = (A > B) ? A : B ;

    if( diff < largest * 1E-9)
    {
        return true ;
    }
    return false ;
}


/* Check that a <= b */
bool lteq
(
    VSPFLOAT a,
    VSPFLOAT b
)
{
    if(a < b)
    {
        return true ;
    }
    else if ( equal(a, b) )
    {
        return true ;
    }
    else
    {
        return false ;
    }
}

// Return true if a < b, and false otherwise
bool Is_Less_Than(VSPFLOAT a, VSPFLOAT b, VSPFLOAT tol)
{
    return a < b - tol ;
}

// Return true if a > b, and false otherwise
bool Is_Greater_Than(VSPFLOAT a, VSPFLOAT b, VSPFLOAT tol)
{
    return a > b + tol ;
}

VSPFLOAT Azi
(
    VSPINT i,
    Az *Az
)
{
    VSPFLOAT Azi = Az->int_part[i] ;
    VSPINT ind_frac_z = Az->ind_frac_z ;

    if(ind_frac_z != EMPTY)
    {
        if( Az->col_frac_z[i] == Az->one )
        Azi += Az->val_frac_z ;
    }
    return Azi ;
}

/* Compute c[i] - gamma * Az[i] */
VSPFLOAT Gradi(VSPINT i, VSPINT *c, VSPFLOAT gamma, Az *Az)
{
    return c[i] - gamma * Azi(i, Az) ;
}




// Round down the fractional component, if there is one
void State_Round_Down
(
    State         *Z,
    Problem *problem
)
{
    VSPINT k = Z->point->ind_frac ;
    if(k == EMPTY) return ;

    VSPINT *c = problem->graph->c ;

    // Reduce the cost
    Z->weight = ((VSPFLOAT) Z->weight) - c[k] * Z->point->val_frac ;

    // Eliminate the fractional entry from z
    Z->point->ind_frac = EMPTY ;
    Z->point->val_frac = 0 ;

    // Eliminte the fractional entry from Az
    Z->Az->ind_frac_z = EMPTY ;
    Z->Az->val_frac_z = 0 ;

    #ifdef DEBUG
        check_state(Z, problem) ;
    #endif

}

void Round_Down_Fractions
(
    State           *Z,
    State           *W,
    Problem   *problem
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        check_states(Z, W, problem) ;
    #endif

    State_Round_Down(Z, problem) ;
    State_Round_Down(W, problem) ;

    #ifdef DEBUG
        check_states(Z, W, problem) ;
    #endif
}


void Az_pm_col_of_A
(
    VSPINT             *Az, // vector to which we add/subtract
    VSPINT               i, // column to add/subtract
    VSPINT          operation, // 1 if add, -1 if subtract
    Problem   *problem
)
{
    VSPINT *Ai = problem->graph->Ai ;
    VSPLONG *Ap = problem->graph->Ap ;
    VSPLONG const end = Ap [i+1] ;
    for(VSPLONG p = Ap[i] ; p < end ; p++)
    {
        VSPINT const j = Ai [p] ;
        Az[j] += operation;
    }
    // don't forget the identity matrix
    Az[i] += operation;
}

/* Update Az struct, given z as a point struct */
void Compute_Az
(
    Point              *z,
    Az                *Az,
    VSPINT                *Ai,
    VSPLONG               *Ap,
    VSPINT                  n
)
{
    SET_FUNCTION_COUNTER ;

    // Compute integer part
    mat_vec_mult(z->int_part, Az->int_part, Ai, Ap, n) ;
    Az->ind_frac_z = z->ind_frac ;

    // Compute fractional part
    Az->ind_frac_z = z->ind_frac ;
    Az->val_frac_z = z->val_frac ;

    if (z->ind_frac != EMPTY)
    {
        Az->one++ ;
        extract_matrix_col(Az->col_frac_z, z->ind_frac, Az->one, Ai, Ap) ;
    }
}

void Error
( 
    char* msg,
    char* file_name,
    VSPINT line_num,
    unsigned long long count
)
{
    #ifndef DEBUG
        printf("\n\nWARNING! DOING DEBUGGING CHECKS WHILE NOT IN DEBUG MODE") ;
    #endif

    #ifdef DEBUG
        static long long unsigned error_count = 0 ;
        error_count++ ;


        FILE* file = fopen ( "output/Error_log.txt", "a" ) ;
        if(error_count <= 10000)
        {
            fprintf(file, "\n\n ERROR: %s", msg) ;
            fprintf(file, "\n File name   : %s", file_name) ;
            fprintf(file, "\n Line number : %lld", (long long) line_num) ;
            fprintf(file, "\n Count : %llu", count) ;
            fclose(file) ;
        }
        else
        {
            fprintf(file, "\n\n: Error log too large, terminating program.\n") ;
            fclose(file) ;
            printf("\n\nERROR: Error log too large, terminating program.\n\n") ;
            exit(1) ;
        }
    #endif
    return ;
}


// Given the Ai, Ap arrays of a binary matrix, extract its kth Ak column into an
// array, with the argument "one" representing 1's
void extract_matrix_col
(
    VSPINT           *Ak, /* Array to store the column */
    VSPINT             k, /* Index of column to extract */
    VSPINT           one, /* Representative of 1's */
    VSPINT           *Ai,
    VSPLONG          *Ap
) {
    VSPINT  j ;
    VSPLONG p, end ;

    end = Ap [ k + 1 ] ;
    for ( p = Ap [k]; p < end; p++ ) {
        j = Ai [p] ;
        Ak [j] = one ;
    }

    /* Take identity into account */
    Ak [k] = one ;

    return ;
}


long double compute_elapsed_time(struct timespec start, struct timespec end)
{
    return    CAST( (end.tv_sec - start.tv_sec), long double) + CAST( (end.tv_nsec - start.tv_nsec), long double) / CAST(1e9, long double) ;
}


VSPFLOAT* Az_struct_to_array (Az *Az) {
    VSPINT    i ;
    VSPFLOAT  frac_entry ;
    VSPINT    ind_frac_z = Az->ind_frac_z ;
    VSPINT   *Az_int = Az->int_part ;
    VSPINT    n      = Az->n ;

    VSPFLOAT *Az_float = malloc(n * sizeof(*Az_float)  ) ;
    VSPFLOAT  val_frac_z = Az->val_frac_z ;


    for (i = 0 ; i < n ; i++) Az_float[i] = (VSPFLOAT) Az_int[i] ;

    // Sometimes col_frac_z is null when we make states in debugging
    if(Az->col_frac_z != NULL)
    {
        for(i = 0 ; i < n ; i++)
        {
            frac_entry = (ind_frac_z != EMPTY) * val_frac_z ;
            Az_float[i] += frac_entry * ( Az->col_frac_z[i] == Az->one ) ;
        }
    }

    return Az_float ;
}


void array_to_point(
    VSPFLOAT  *z_float,
    Point       *z
) {
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
        bool frac_entry = false ;
    #endif
    VSPINT i ;
    VSPINT n = z->n ;
    VSPINT *z_int = z->int_part ;

    z->ind_frac = EMPTY ;
    for(i = 0 ; i < n ; i++) {
        if( fabs(z_float[i]) < 1E-6 )            z_int[i] = 0 ;
        else if ( fabs(z_float[i] - 1) < 1E-6 )  z_int[i] = 1 ;
        else {
            #ifdef DEBUG
                if (frac_entry) {
                    Error("Multiple frac entries", DBG_INFO) ;
                }
                frac_entry = true ;
            #endif
            z_int[i]    = 0 ;
            z->ind_frac = i ;
            z->val_frac = z_float[i] ;
        }
    }
    return ;
}

VSPFLOAT* Point_to_array(Point *z) {
    VSPINT i, k ;
    VSPINT n = z->n ;

    /* Allocate space */
    VSPFLOAT *z_arr = malloc(n * sizeof(*z_arr)) ;

    /* Copy the integer part */
    for(i = 0 ; i < n ; i++) z_arr[i] = z->int_part[i] ;

    k = z->ind_frac ;
    if(k != EMPTY) {
        z_arr[k] = z->val_frac ;
    }

    return z_arr ;
}



// Create a deep copy. Assume that all memory is already allocated
void copy_state(State *from, State *to)
{
    VSPINT  i ;
    VSPINT  n = from->point->n ;
    VSPINT *from_z_int  = from->point->int_part ;
    VSPINT *to_z_int    = to  ->point->int_part ;
    VSPINT *from_Az_int = from->Az->int_part ;
    VSPINT *to_Az_int   =   to->Az->int_part ;

    /* Copy n */
    /* TODO speed: The way I'm using copy_state, copying the n should not
       be necessary. I am only doing it here to be safe                       */
    to->point->n = from->point->n ;
    to->Az   ->n = from->   Az->n ;

    /* Copy the weight */
    to->weight = from->weight ;

    /* If from has a fractional entry */
    if(from->point->ind_frac != EMPTY) {
        VSPINT *from_col_A  = from->Az->col_frac_z ;
        VSPINT *to_col_A    = to->Az->col_frac_z ;
        VSPINT  from_one    = from->Az->one ;

        to->Az->one++ ;
        VSPINT  to_one      = to->Az->one ;

        /* Copy all the arrays in the state struct */
        for(i = 0 ; i < n ; i++) {
            to_z_int [i] = from_z_int[i] ;
            to_Az_int[i] = from_Az_int[i] ;
            /* TODO speed: perhaps branchless code is not the best here.      */
            to_col_A [i] = to_one * (from_col_A[i] == from_one) ;
        }

        /* Copy the rest of the info about the fractional entry */
        to->point->ind_frac = from->point->ind_frac ;
        to->point->val_frac = from->point->val_frac ;
        to->Az->ind_frac_z  = from->Az->ind_frac_z ;
        to->Az->val_frac_z  = from->Az->val_frac_z ;
    }
    /* If from does not have a fractional entry */
    else {
        /* Copy all the arrays in the state struct */
        for(i = 0 ; i < n ; i++) {
            to_z_int [i] = from_z_int[i] ;
            to_Az_int[i] = from_Az_int[i] ;

            /* Copy the rest of the info about the fractional entry */
            to->point->ind_frac = from->point->ind_frac ;
            to->point->val_frac = from->point->val_frac ;
            to->Az->ind_frac_z  = from->Az->ind_frac_z ;
            to->Az->val_frac_z  = from->Az->val_frac_z ;
        }
    }
    return ;
}

// Allocate + copy
State *Clone_State(State *src)
{
    State *dst = Construct_State(src->point->n);
    copy_state(src, dst);
    return dst;
}



/* Use when you want to sort an array in increasing order using qsort from
   stdlib                                                                     */
int comparison_increasing(const void *a, const void *b) {
    if ( *(VSPINT*) a < *(VSPINT*) b ) {
        return -1 ;
    }
    else if ( *(VSPINT*) a > *(VSPINT*) b ) {
        return 1 ;
    }
    else return 0 ;
}

/* Dump the point in the file "output/title" */
void DumpPoint(Point *z, char *outName)
{
    /* Open a file in output */
    char fileName[300] = {'\0'} ;
    strcpy(fileName, "output/") ;
    strcat(fileName, outName) ;
    FILE *file = fopen(fileName, "w") ;


    /* Write each entry on a new line */
    for(VSPINT i = 0 ; i < z->n ; i++ ) 
    {
        if(i == z->ind_frac)
        {
            fprintf(file, "%f\n", z->val_frac) ; 
        }
        else
        {
           fprintf(file, "%lld\n", (long long) z->int_part[i]) ; 
        }
    }
    fclose(file) ;
    return ;
}

// Dump the point in the file "output/title
// File name should not be super long
// The output looks like:
// x,y
// 1,0
// 1,1
// 0,1
void DumpPoints(Point *x, Point *y, char *outName)
{
    // Open a file in output
    char fileName[320] = {'\0'} ;
    strcpy(fileName, "output/") ;
    strcat(fileName, outName) ;
    FILE *file = fopen(fileName, "w") ;


    // x
    fprintf(file, "x,y") ;
    for(VSPINT i = 0 ; i < x->n ; i++ ) 
    {
        // x
        fprintf(file, "\n") ;
        if(i == x->ind_frac)
        {
            fprintf(file, "%f,", x->val_frac) ; 
        }
        else
        {
           fprintf(file, "%lld,", (long long) x->int_part[i]) ; 
        }

        // y
        if(i == y->ind_frac)
        {
            fprintf(file, "%f", y->val_frac) ; 
        }
        else
        {
           fprintf(file, "%lld", (long long) y->int_part[i]) ; 
        }
    }

    fclose(file) ;
    return ;
}

void DumpIntArray(VSPINT *array, VSPINT len, char *outName)
{
    /* Open a file in output */
    char fileName[300] = {'\0'} ;
    strcpy(fileName, "output/") ;
    strcat(fileName, outName) ;
    FILE *file = fopen(fileName, "w") ;

    for(VSPINT i = 0 ; i < len ; i++)
    {
        fprintf(file, "%lld  ", (long long) array[i]) ;
    }
    fclose(file) ;
    return ;
}


/* Set frac entry, but do NOT update the weight */
void SetFracEntry
(
    VSPINT         indFrac,
    VSPFLOAT       valFrac,
    State           *Z,
    Problem   *problem
)
{
    #ifdef DEBUG
        check_state(Z, problem) ;
    #endif
    Z->point->ind_frac = indFrac ;
    Z->point->val_frac = valFrac ;

    Z->Az->ind_frac_z = indFrac ;
    Z->Az->val_frac_z = valFrac ;

    Z->Az->one++ ;
    extract_matrix_col
    (
        Z->Az->col_frac_z,
        indFrac,
        Z->Az->one,
        problem->graph->Ai,
        problem->graph->Ap
    ) ;

    #ifdef DEBUG
        check_state(Z, problem) ;
    #endif
    return ;
}

VSPFLOAT ComputeFracEntry
(
    VSPINT   u,
    VSPINT wgt,
    VSPINT  ci
)
{
    return (u - wgt) / ( (VSPFLOAT) ci ) ;
}


VSPINT Count_Ones
(
    VSPINT *x, // assumed binary
    VSPINT  n  // size of x
)
{
    VSPINT ret = 0 ;
    for(VSPINT i = 0 ; i < n ; i++) ret += x[i] ;
    return ret ;
}

// Returns a permutation of the integers 0, 1, 2, ..., n-1 in arr
// ASSUMED: No assumptions on arr, it can be garbage
void Get_Permutation(VSPINT *arr, VSPINT n)
{
    for (VSPINT i = 0; i < n; i++) arr[i] = i;

    for (VSPINT i = n - 1; i > 0; i--)
    {
        VSPINT j = rand() % (i + 1);
        VSPINT temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}


void Swap_Mem(void *a, void *b, size_t size)
{
    void *temp = malloc(size);  

    memcpy(temp, a, size);
    memcpy(a, b, size);
    memcpy(b, temp, size);

    free(temp); // Free the temporary buffer
}

// Return a pointer to a state containing a deep copy of z, Az, and wgtz
void Save_State_Copy
(
    VSPINT *z, VSPINT *Az, VSPINT wgtz, VSPINT n, State **stateHome
)
{
    State *ret = Construct_State(n) ;
    copy_array_int(z , ret->point->int_part, n) ;
    copy_array_int(Az, ret->Az   ->int_part, n) ;
    ret->weight = wgtz ;

    *stateHome = ret ;
    return ;
}

void Check_Array_Is_Binary
(
    VSPINT *array,
    VSPINT      n
)
{
    static long long unsigned count = 0 ;
    count++ ;
    VSP_ASSERT(n > 0) ;

    VSPINT i ;
    for(i = 0 ; i < n ; i++)
    {
        if(array[i] != 0 && array[i] != 1)
        {
            ERR("Non-binary entry") ;
            break ;
        }
    }
}


// base ^ expopnent where both base and exponent are
// 0^0 is 1
VSPINT pow_ints(VSPINT base, VSPINT exponent)
{
    VSPINT ret = 1 ;
    for(VSPINT i = 0 ; i < exponent ; i++)
    {
        ret *= base ;
    }
    return ret ;
}


// (nz/2) / n(n-1)/2 = nz/n(n-1)
double Graph_Density(VSPINT n, VSPLONG nz)
{
    double n_dbl  = CAST(n, double) ;
    double nz_dbl = CAST(nz, double) ;
    return nz_dbl / ( n_dbl * (n_dbl - 1) ) ;
}

VSPINT Get_Max_Degree(VSPLONG *Ap, VSPINT n)
{
    VSPINT i, deg, max_deg = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        deg = Ap[i+1] - Ap[i] ;
        if(deg > max_deg) max_deg = deg ;
    }
    return max_deg ;
}

bool Is_Singleton
(
    VSPINT    i,
    VSPLONG *Ap
)
{
    return (Ap[i+1] - Ap[i] == 0) ;
}


// Write 4n EMPTYs, covering x, y, Ax, Ay
static void Write_EMPTYs_When_NULL
(
    FILE    *fp,
    int32_t   n
)
{
    int32_t empty = EMPTY ;
    for(int32_t i = 0 ; i < 4*n ; i++)
    {
        fwrite(&empty, sizeof(int32_t), 1, fp) ;
    }
}

// problem.hager6 n nz Ai Ap c x y (integer parts) Ax Ay (integer parts)
// problem.txt Contains u, and any fractional entries of X and Y
void Dump_Problem
(
    const Problem *P,
    const State   *X,
    const State   *Y
)
{
    int32_t i_32 ;
    int64_t i_64 ;

    FILE *fp  = fopen("output/problem.hager6", "wb") ;
    FILE *fp2 = fopen("output/problem.txt", "w") ;
    
    // Write n
    int32_t n_32 = (int32_t) P->graph->n ;
    fwrite(&n_32, sizeof(int32_t), 1, fp) ;

    // Write nz
    int64_t nz_64 = (int64_t) P->graph->nz ;
    fwrite(&nz_64, sizeof(int64_t), 1, fp) ;

    //  Write Ai
    VSPINT *Ai = P->graph->Ai ;
    int32_t *Ai_32 = malloc( nz_64 * sizeof(*Ai_32) ) ;
    for(i_64 = 0 ; i_64 < nz_64 ; i_64++) Ai_32[i_64] = (int32_t) Ai[i_64] ;
    fwrite(Ai_32, sizeof(int32_t), nz_64, fp) ;
    free(Ai_32) ;
    Ai_32 = NULL ;

    //  Write Ap
    VSPLONG *Ap = P->graph->Ap ;
    int64_t *Ap_64 = malloc( (n_32 + 1) * sizeof(*Ap_64) ) ;
    for(i_32 = 0 ; i_32 < n_32 + 1 ; i_32++) Ap_64[i_32] = (int64_t) Ap[i_32] ;
    fwrite(Ap_64, sizeof(int64_t), n_32 + 1, fp) ;
    free(Ap_64) ;
    Ap_64 = NULL ;

    // Write c
    VSPINT *c = P->graph->c ;
    int32_t *c_32 = malloc( n_32 * sizeof(*c_32) ) ;
    for(i_32 = 0 ; i_32 < n_32 ; i_32++) c_32[i_32] = (int32_t) c[i_32] ;
    fwrite(c_32, sizeof(int32_t), n_32, fp) ;
    free(c_32) ;
    c_32 = NULL ;

    if(X != NULL && Y != NULL)
    {
        // Write x
        VSPINT *x = X->point->int_part ;
        int32_t *x_32 = malloc( n_32 * sizeof(*x_32) ) ;
        for(i_32 = 0 ; i_32 < n_32 ; i_32++) x_32[i_32] = (int32_t) x[i_32] ;
        fwrite(x_32, sizeof(int32_t), n_32, fp) ;
        free(x_32) ;
        x_32 = NULL ;

        // Write y
        VSPINT *y = Y->point->int_part ;
        int32_t *y_32 = malloc( n_32 * sizeof(*y_32) ) ;
        for(i_32 = 0 ; i_32 < n_32 ; i_32++) y_32[i_32] = (int32_t) y[i_32] ;
        fwrite(y_32, sizeof(int32_t), n_32, fp) ;
        free(y_32) ;
        y_32 = NULL ;

        // Write Ax
        VSPINT *Ax = X->Az->int_part ;
        int32_t *Ax_32 = malloc( n_32 * sizeof(*Ax_32) ) ;
        for(i_32 = 0 ; i_32 < n_32 ; i_32++) Ax_32[i_32] = (int32_t) Ax[i_32] ;
        fwrite(Ax_32, sizeof(int32_t), n_32, fp) ;
        free(Ax_32) ;
        Ax_32 = NULL ;

        // Write Ay
        VSPINT *Ay = Y->Az->int_part ;
        int32_t *Ay_32 = malloc( n_32 * sizeof(*Ay_32) ) ;
        for(i_32 = 0 ; i_32 < n_32 ; i_32++) Ay_32[i_32] = (int32_t) Ay[i_32] ;
        fwrite(Ay_32, sizeof(int32_t), n_32, fp) ;
        free(Ay_32) ;
        Ay_32 = NULL ;

        // Write fractional entries of x and y
        fprintf(fp2, "Ind frac x: %lld\n", (long long) X->point->ind_frac) ;
        fprintf(fp2, "Val frac x: %f\n"  ,  X->point->val_frac) ;
        fprintf(fp2, "\nInd frac y: %lld\n", (long long) Y->point->ind_frac) ;
        fprintf(fp2, "Val frac y: %f\n"  ,  Y->point->val_frac) ;
    }
    else
    {
        Write_EMPTYs_When_NULL(fp, n_32) ;
    }
    fclose(fp) ;
    fp = NULL ;

    // Write u and level
    fprintf(fp2, "level     : %lld\n", (long long) P->level) ;
    fprintf(fp2, "u         : %lld\n", (long long) P->u) ;
    fclose(fp2) ;
    return ;
}

VSPINT State_Get_Wgt_Int_Part
(
    State *Z,
    VSPINT    *c
)
{
    if(Z->point->ind_frac == EMPTY) return Z->weight ;

    VSPFLOAT frac_wgt = c[Z->point->ind_frac] * Z->point->val_frac ;
    return (VSPINT) (Z->weight - frac_wgt + 1e-3) ;
}

// Assume: array is not empty
VSPINT Get_Rand_Array_Entry
(
    VSPINT *arr,
    VSPINT    n  // size of arr
)
{
    return arr[ rand() % n ] ;
}


/* ************* Explicit Matrix ******************************************* **/
VSPINT** construct_explicit_matrix
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
        Error("Memory allocation failed", DBG_INFO) ;
        exit(1) ;
    }

    /* Allocate space for each row */
    for (i = 0 ; i < n ; i++) {
        A [i] = malloc( n * sizeof(VSPINT) ) ;
        if ( A [i] == NULL ) {
            Error("Memory allocation failed", DBG_INFO) ;
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
            Error("Adjacency matrix not symmetric", DBG_INFO) ;
        }
        else
        {
            Error("Edge weights matrix not symmetric", DBG_INFO) ;
        }
    }

    return A ;
}

void destruct_explicit_matrix( VSPINT** A, VSPINT n ) {
    VSPINT i ;

    for ( i = 0 ; i < n ; i++ ) free (A [i]) ;
    free(A) ;

    return ;
}




void DumpNeighs(VSPINT i, VSPINT *Ai, VSPLONG *Ap)
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

void ResetDumpFile()
{
    FILE *file = fopen("output/dump.txt", "w") ;
    fclose(file) ;
    return ;
}

void check_matrix_is_symmetric(
    VSPINT           *Ai,
    VSPLONG *Ap,
    VSPINT           *Ae,
    VSPLONG nz,
    VSPINT             n,
    Options* options
) ;

/* Check each state individually, as well as the cross term */
void check_states(State *X, State *Y, Problem *problem)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    check_state(X, problem) ;
    check_state(Y, problem) ;

    /* Check the cross term */
    VSPFLOAT *x_f = Point_to_array(X->point) ;
    VSPFLOAT *y_f = Point_to_array(Y->point) ;
    VSPFLOAT *Ax_f = Az_struct_to_array(X->Az) ;
    VSPFLOAT *Ay_f = Az_struct_to_array(Y->Az) ;

    VSPFLOAT xAy_1 = 0 ;
    VSPFLOAT xAy_2 = 0 ;

    for(VSPINT i = 0 ; i < problem->graph->n ; i++)
    {
        xAy_1 += x_f[i] * Ay_f[i] ;
        xAy_2 += y_f[i] * Ax_f[i] ;
    }

    if( equal(xAy_1, xAy_2) == false )
    {
        Error("Cross terms don't agree", DBG_INFO) ;
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
    Problem *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;
    WARN_IF_NOT_DEBUG_MODE ;
    VSPINT n = problem->graph->n ;
    /*------------------------------------------------------------------------*/
    // Load the data into states and use check_states
    /*------------------------------------------------------------------------*/
    Point xPoint, yPoint ;

    xPoint.int_part = x ;
    xPoint.ind_frac = EMPTY ;
    xPoint.val_frac = 0 ;
    xPoint.n = n ;
    
    yPoint.int_part = y ;
    yPoint.ind_frac = EMPTY ;
    yPoint.val_frac = 0 ;
    yPoint.n = n ;
    
    Az Axstruct, Aystruct ;

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

    State X, Y ;

    X.point  = &xPoint ;
    X.Az     = &Axstruct ;
    X.weight = wgtx ;

    Y.point  = &yPoint ;
    Y.Az     = &Aystruct ;
    Y.weight = wgty ;

    check_states(&X, &Y, problem) ;


    if(wgtx < 0 || wgty < 0)
    {
        ERR("Negative weights") ;
    }

    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;
    return ;
}

void check_max_of_array
(
    VSPINT        *c,
    VSPINT         n, /* Length of c */
    VSPINT max_given /* Purported maximum value of c */
) ;


/* Checks:
    - The point is feasible
    - Az is correct
    - The weight is correct   */
void check_state(State *Z, Problem *problem)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    /* Check that the n's agree */
    VSPINT nFromZ       = Z->point->n ;
    VSPINT nFromAz      = Z->Az->n ;
    VSPINT nFromProblem = problem->graph->n ;
    if(nFromZ != nFromProblem)
    {
        Error("n is incorrect", DBG_INFO) ;
    }

    if(nFromAz != nFromProblem)
    {
        Error("n is incorrect", DBG_INFO) ;
    }

    /* Check the point */
    check_point(Z->point) ;

    /* Check that Az is correct */
    check_Az(Z->point, Z->Az, problem) ;

    /* Check that the weight is correct, and z satisfies constraint */
    check_weight(Z, problem->graph->c) ;
    if(Z->weight > problem->u)
    {
        Error("Upper balance violated ", DBG_INFO);
    }

    return ;
}



void check_weight(State *Z, VSPINT *c) {

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

    if( !equal(weight_check_int + weight_check_frac, Z->weight) )
    {
        Error("Incorrect weight", DBG_INFO) ;
    }
}

void check_Az
(
    Point         *z,
    Az           *az,
    Problem *problem
)
{

    static long long unsigned count = 0 ;
    count ++ ;

    VSPINT i ;
    VSPINT n       = problem->graph->n ;
    VSPINT *Az_int = az->int_part ;

    /* Compute Az from scratch */
    Az *Az_dbg = Construct_Az(n) ;
    Compute_Az(z, Az_dbg, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;

    /* Check that each entry in the integer part is in [0, n] */
    for(i = 0 ; i < n ; i++) {
        if(Az_int[i] < 0 || Az_int[i] > n) {
            Error("Problem in Az", DBG_INFO) ;
            break ;
        }
    }

    /* Check the integer parts are equal */
    check_arrays_are_equal(az->int_part, Az_dbg->int_part, n) ;

    /* Check that the fractional part is equal */
    if(az->ind_frac_z != Az_dbg->ind_frac_z) {
        Error("Error in Az", DBG_INFO) ;
    }

    if(az->ind_frac_z != EMPTY) {

        /* Check that the integer part of the fractional entry is < n*/
        if(Az_int[az->ind_frac_z] > n - 1 || Az_int[az->ind_frac_z] < 0) {
            Error("Problem with fractional index in Az", DBG_INFO) ;
        }

        /* Check that the fractional columns are the same */
        bool    same              = true ;
        bool    one_is_ok         = true ;
        VSPINT *col_frac_z        = az->col_frac_z ;
        VSPINT *col_frac_z_dbg    = Az_dbg->col_frac_z ;
        VSPINT  one               = az->one ;
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
            Error("Fractional columns mismatch", DBG_INFO) ;
        }
        if(one_is_ok == false) {
            Error("Problem with fractional column", DBG_INFO) ;
        }

        /* Check that the fractional column has a one in index ind_frac_z */
        if(col_frac_z[az->ind_frac_z] != one) {
            Error("Matrix column wrong", DBG_INFO) ;
        }
    }

    /* Free memory */
    Destroy_Az(Az_dbg) ;

    return ;
}



void check_arrays_are_equal
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
        Error("Arrays not equal", DBG_INFO) ;
    }
    return ;
}

void Check_Arrays_Are_Equal_VSPLONG
(
    const VSPLONG *const u, /* Check that vector u equals vector v */
    const VSPLONG *const v,
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
        Error("Arrays not equal", DBG_INFO) ;
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
VSPINT check_matrix
(
    VSPINT* Ai,
    VSPLONG* Ap,
    VSPINT* Ae,
    VSPLONG nz,
    VSPINT n,
    Options* options
) {
    static long long unsigned count = 0 ;
    count += 1;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, j, k ;
    VSPLONG end, p ;
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
                Error(" Error in edge weights", DBG_INFO) ;
            }
            if ( Ae [p] > 1E10 ) {
                Error("Warning: probable error in edge weights", DBG_INFO) ;
            }

        }

        /* Check that the diagonal is 0 */
        if (col [k] != 0) {
            diag_nonzero = TRUE ;
        }
    }

    /* Check that the matrix is symmetric */
    check_matrix_is_symmetric(Ai, Ap, Ae, nz, n, options) ;

    /* Print errors, if there are any */
    if ( Ap_not_increasing ) {
        Error("Ap not increasing", DBG_INFO) ;
    }

    if ( repeated_entry ) {
        Error("Repeated neighbor in Ai", DBG_INFO) ;
        error = TRUE ;
    }

    if ( diag_nonzero ) {
        Error("Nonzero diagonal", DBG_INFO) ;
        error = TRUE ;
    }

    if ( invalid_neighbor ) {
        Error("Invalid neighbor in graph", DBG_INFO) ;
        error = TRUE ;
    }

    if ( last_entry_of_Ap_incorrect ) {
        Error("Last entry of Ap incorrect", DBG_INFO) ;
    }

    if ( invalid_entry_in_Ap ) {
        Error("Invalid entry in Ap", DBG_INFO) ;
    }

    /* Free memory */
    free(col) ;

    return error ;
}


void check_matrix_is_symmetric
(
    VSPINT           *Ai,
    VSPLONG *Ap,
    VSPINT           *Ae,
    VSPLONG nz,
    VSPINT             n,
    Options* options
) {
    /* Get the transpose */
    Problem* T = transpose(Ai, Ap, Ae, nz, n, options) ;
    VSPINT *Ti = T->graph->Ai, *Te = T->graph->Ae ;
    VSPLONG *Tp = T->graph->Ap ;

    /* Check that the number of 1's in each column is the same in A and T */
    Check_Arrays_Are_Equal_VSPLONG(Ap, Tp, n+1) ;

    /* Check that the row index of each 1 in each column is the same in A, T*/
    for(VSPINT j = 0 ; j < n ; j++) {
        check_permutation(&Ai[Ap[j]], &Ti[Tp[j]], Ap[j+1] - Ap[j]) ;
    }

    if (Ae != NULL) check_permutation(Ae, Te, nz) ;

    /* Free memory */
    Destroy_Problem(T) ;
}


/* ============================================================================
   ================ check_entries_are_between_bounds ======================
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

void check_entries_are_between_bounds(
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
        Error("Array entry out of bounds", DBG_INFO) ;
    }

}

/* ============================================================================
   ================ check_max_of_array ====================================
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

void check_max_of_array(
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
        Error("Max of array not correct", DBG_INFO) ;
    }
}

/* ============================================================================
   ================ check_problem =========================================
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
 a) A Problem struct

 OUTPUTS:
 a) If the matrix passes all the checks, return TRUE (i.e. integer 1).
    Otherwise, return FALSE (integer 0) and print errors to the error log.
*/

VSPINT check_problem
(
    Problem* problem,
    Options* options
)
{
    static long long unsigned count = 0 ;
    count += 1 ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i, deg ;

    /* Unpack problem */
    VSPINT* Ai = problem->graph->Ai ;
    VSPLONG* Ap = problem->graph->Ap ;
    VSPINT* Ae = problem->graph->Ae ;
    VSPINT* c = problem->graph->c ;

    VSPINT alpha = problem->graph->min_vtx_cost ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPLONG nz    = problem->graph->nz ;
    VSPINT n     = problem->graph->n ;
    VSPINT bal   = problem->u ;

    VSPINT cost_is_wrong     = FALSE ;
    VSPINT weight_is_wrong   = FALSE ;
    VSPINT is_there_an_error = FALSE ;


    /* Do all the matrix checks */
    check_matrix(Ai, Ap, Ae, nz, n, options) ;

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

    if (cost_is_wrong) {
        Error("Cost is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (weight_is_wrong) {
        Error("Weight is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (locally_computed_alpha != alpha) {
        Error("Gamma is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (locally_computed_gamma != gamma) {
        Error("Gamma is wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (bal < 1) {
        Error("Balance constraints are wrong", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    if (sum_of_all_ci != problem->graph->total_vtx_cost) {
        Error("Problem with sum of all ci", DBG_INFO) ;
        is_there_an_error = TRUE ;
    }

    return is_there_an_error ;

}



/* Check whether arr_2 is a permutation of arr_1 */
void check_permutation (
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
    qsort(arr_1_copy, n, sizeof(*arr_1), comparison_increasing) ;
    qsort(arr_2_copy, n, sizeof(*arr_2), comparison_increasing) ;

    VSPINT arrays_are_equal = TRUE ;
    for (i = 0 ; i < n ; i++) {
        if (arr_1_copy [i] != arr_2_copy [i]) {
            arrays_are_equal = FALSE ;
        }
    }

    if (arrays_are_equal == FALSE) {
        Error("Array not a permutation", DBG_INFO ) ;
    }

    /* Free memory */
    free(arr_1_copy) ;
    free(arr_2_copy) ;
    return ;
}


/* Check as much as possible when only the Az struct is available             */
void light_check_Az(
    Az *Az
) {

    static long long unsigned count = 0 ;
    count++ ;

    VSPINT  n      = Az->n ;
    VSPINT *Az_int = Az->int_part ;

    /* Check that the integer part is between 0 and n */
    for(VSPINT i = 0 ; i < n ; i++) {
        if(Az_int[i] < 0 || Az_int[i] > n) {
            Error("Error in Az", DBG_INFO) ;
        }
    }


    /* Check that the fractional column has a one in index ind_frac_z */
    if(Az->ind_frac_z != EMPTY)
        {
        if(Az->col_frac_z[Az->ind_frac_z] != Az->one) {
            Error("Matrix column wrong", DBG_INFO) ;
        }
    }
}


void check_point(Point *z) {

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
    if( binary == false ) Error("Not binary", DBG_INFO) ;


    /* Check that the integer part of the fractional entry is 0 */
    VSPINT ind_frac = z->ind_frac ;
    if(ind_frac != EMPTY && z_int[ind_frac] != 0) {
        Error("Problem with integer part", DBG_INFO) ;
    }

    /* Check that the fractional index is empty or in the valid range */

    if( ind_frac != EMPTY && !(0 <= ind_frac && ind_frac <= n-1) ) {
        Error("Fractional index out of bounds", DBG_INFO) ;
    }

    /* Check that the fractional value is in [0, 1] */
    VSPFLOAT val_frac = z->val_frac ;
    if(val_frac < 0 || val_frac > 1) {
        Error("Issue with fractional value", DBG_INFO) ;
    }

    /* If frac index exists, check that the frac value is in (0, 1)*/
    if(ind_frac != EMPTY) {
        if(val_frac <= 0 || val_frac >= 1 ) {
            Error("Issue with fractional value", DBG_INFO) ;
        }
    }
    return ;
}



void GetSepCost
(
    State          *X,
    State          *Y,
    Problem  *problem,
    VSPINT       *sepCost   /* Output*/
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif

    VSPINT  n = X->point->n ;
    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;
    VSPINT *c = problem->graph->c ;

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
void xAy
(
    State         *X,
    State         *Y,
    Problem *problem,
    VSPFLOAT        *xAy
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif
    VSPINT i ;
    VSPINT n = problem->graph->n ;
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
void compute_gradient_v1
(
    VSPFLOAT *grad_z,
    Az       *Aw,
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
void mat_vec_mult
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







/* Compute the weight of a point. This function assumes that the true weight
   is an integer.                                                             */
void compute_weight(Point *Z, VSPINT *weight_z, VSPINT *c) {
    #ifdef DEBUG
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
        #ifdef DEBUG
            if( fabs(weight_frac_part - round(weight_frac_part)) > ZERO ) {
                Error("Error with weight", DBG_INFO) ;
            }
        #endif
    }
    *weight_z = weight ;
}

Problem* transpose(
    VSPINT           *Ai, /* Row indices of non-zeros  */
    VSPLONG *Ap, /* Column pointers           */
    VSPINT           *Ae, /* Entries                   */
    VSPLONG nz, /* Number of rows/columns    */
    VSPINT             n, /* Number of nonzeros        */
    Options* options
) {
 
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    #ifdef DEBUG
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
    Problem* T = Construct_Problem(n, nz, NULL, NULL, NULL, options) ;
    VSPINT *Ti = T->graph->Ai , *Te = T->graph->Ae ;
    VSPLONG *Tp = T->graph->Ap ;

    /* Count the number of nonzeros in each row */
    for (VSPLONG ii = 0 ; ii < nz ; ii++) {
        ind_row = Ai [ii] ;
        num_nonzeros_in_row [ ind_row ] += 1 ;
    }


    /* Set Tp[0] = 0 and Tp[i] = number of nonzeros in rows 0, 1, ..., i-1 */
    /* Then, set num_nonzeros_in_row [i] = Tp [i] for i = 0, ..., n-1 */
    partial_sums_of_array(num_nonzeros_in_row, Tp, n) ;

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

    #ifdef DEBUG

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
            VSPINT** A_explicit = construct_explicit_matrix(NULL, Ai, nz,
                                                                      Ap, n  ) ;
            VSPINT** T_explicit = construct_explicit_matrix(NULL, Ti, nz,
                                                                      Tp, n  ) ;

            VSPINT** A_vals_explicit = NULL ;
            VSPINT** T_vals_explicit = NULL ;

            if (Ae != NULL) {
                A_vals_explicit = construct_explicit_matrix(
                                                            Ae, Ai, nz, Ap, n) ;

                T_vals_explicit = construct_explicit_matrix(
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
            destruct_explicit_matrix(A_explicit, n) ;
            destruct_explicit_matrix(T_explicit, n) ;

            if (Ae != NULL) {
                destruct_explicit_matrix(A_vals_explicit, n) ;
                destruct_explicit_matrix(T_vals_explicit, n) ;
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
VSPLONG partial_sums_of_array(
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
VSPFLOAT compute_obj(
    State         *X,
    State         *Y,
    Problem *problem,
    VSPFLOAT       gamma
) {
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        /* Check input */
        check_states(X, Y, problem) ;
    #endif

    /* Local vars */
    VSPINT i ;
    VSPINT n = problem->graph->n ;
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

    #ifdef DEBUG
        if(n < 200) {
            VSPFLOAT obj_val_dbg ;
            obj_val_dbg = compute_objective_0(X, Y, problem, gamma) ;
            if( fabs(obj_val - obj_val_dbg) > ZERO )
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
VSPFLOAT compute_objective_0(
    State *X,
    State *Y,
    Problem *problem,
    VSPFLOAT gamma
) {
    /* Local vars */
    VSPINT i, j ;
    VSPLONG start, end, p ;
    VSPFLOAT  obj_val = 0 ;
    VSPINT    n       = problem->graph->n ;
    VSPINT   *Ai      = problem->graph->Ai ;
    VSPLONG *Ap      = problem->graph->Ap ;
    VSPINT   *c       = problem->graph->c ;


    /* Convert X, Y to arrays of floats */
    VSPFLOAT *x_arr = Point_to_array(X->point) ;
    VSPFLOAT *y_arr = Point_to_array(Y->point) ;

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
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c = problem->graph->c ;

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
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;

    objval = 0 ;


    
    for(i = 0 ; i < n ; i++)
    {
        objval += ( c[i] * (x[i] + y[i]) - gamma * x[i] * Ay[i] ) ;
    }

    #ifdef DEBUG
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
VSPFLOAT compute_LP_objective(
    Point    *z,
    VSPINT       *c,
    VSPFLOAT  gamma,
    Az      *Aw
) {
    #ifdef DEBUG
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
        Awi = Azi(i, Aw) ;
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
void ComputeGradFlt
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
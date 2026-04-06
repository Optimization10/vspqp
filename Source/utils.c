#include "core.h"
#include "shared_api.h"
#include "ds_api.h"
#include <stdint.h>
#include <string.h>
#include <stdarg.h>
#include "flips_common_api.h"



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



void resolve_frac_entry
(
    State *Z,
    State *W,
    Problem *problem
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        check_state(Z, problem) ;
    #endif

    VSPINT k = Z->point->ind_frac ;
    if(k == EMPTY) return ;

    VSPFLOAT val_frac = Z->point->val_frac ;
    VSPINT         *c = problem->graph->c, j ;
    VSPLONG p ;

    /* Check if setting frac entry to 0 would violate lower bound */
    if(Z->weight - val_frac * c[k] < 1 - ZERO)
    {
        /* Round up */
        if(Z->weight + (1 - val_frac) * c[k] <= problem->u)
        {
            VSPLONG *Ap = problem->graph->Ap ;
            VSPINT *Ai = problem->graph->Ai ;

            Z->point->int_part[k] = 1 ;
            Z->weight += lround((1 - val_frac) * c[k]) ;
            Z->point->ind_frac = EMPTY ;

            Z->Az->ind_frac_z = EMPTY ;
            Az_pm_col_of_A(Z->Az->int_part, k, PLUS, problem) ;
            
            // Downflip the neighbors of k in the other variable as long as it
            // doesn't violate the lower balance constraint
            VSPINT *w = W->point->int_part ;
            VSPINT *Aw = W->Az->int_part ;
            for(p = Ap[k] ; p < Ap[k+1] ; p++)
            {
                j = Ai[p] ;
                if(w[j] && c[j] > W->weight)
                {
                    Do_Flip_Basic(w, Aw, j, &(W->weight), c, Ai, Ap) ;
                }
            }
            

            return ;
        }
        else
        {
            // printf("\n\nERROR: Could not deal with fractional component. Problem: %s\n\n", problem->name) ;
        }
    }

    /* If we made it past the if statement, we default to rounding down */
    /* Update the weight */
    Z->weight -= lround( Z->point->val_frac * c[k] ) ;

    /* Update the point struct */
    Z->point->ind_frac = EMPTY ;

    /* Update the Az struct */
    Z->Az->ind_frac_z = EMPTY ;

    #ifdef DEBUG
        check_state(Z, problem) ;
    #endif


    return ;
}

void Force_To_Binary_Point
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

// void DumpGrad(VSPINT *c, VSPFLOAT gamma, Az *Az, char *outName)
// {
//     /* Open a file in output */
//     char fileName[300] = {'\0'} ;
//     strcpy(fileName, "output/") ;
//     strcat(fileName, outName) ;
//     FILE *file = fopen(fileName, "w") ;

//     for(VSPINT i = 0 ; i < Az->n ; i++)
//     {
//         VSPFLOAT Azi = Azi(i, Az) ;
//         fprintf(file, "%f\n", c[i] - gamma * Azi) ;
//     }
//     fclose(file) ;
//     return ;
// }

// void DumpAz(Az *Az, char *outName)
// {
//     /* Open a file in output */
//     char fileName[300] = {'\0'} ;
//     strcpy(fileName, "output/") ;
//     strcat(fileName, outName) ;
//     FILE *file = fopen(fileName, "w") ;

//     for(VSPINT i = 0 ; i < Az->n ; i++)
//     {
//         VSPFLOAT Azi = Azi(i, Az) ;
//         fprintf(file, "%f\n", Azi) ;
//     }
//     fclose(file) ;
//     return ;
// }

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
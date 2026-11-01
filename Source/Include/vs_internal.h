#ifndef _VS_INTERNAL_H_
#define _VS_INTERNAL_H_

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include </usr/include/mkl/mkl.h>
#include </usr/include/mkl/mkl_spblas.h>
#include "vs.h"

/* ========================================================================== */
/* === macros =============================================================== */
/* ========================================================================== */

#define VSINT    VSPINT
#define VSLINT   VSPLONG
#define VSFLOAT  VSPFLOAT

// #define PRI(A) ( (long long) (A) ) // defined in core.h

#define VSFALSE  ( 0)
#define VSTRUE   ( 1)
#define VSEMPTY  (-1)

#define VSONE ((VSFLOAT) 1)
#define VSZERO ((VSFLOAT) 0)
#define VSHALF ((VSFLOAT) .5)
#define RE_SCALAR ((VSFLOAT) 5.85498315243191680e-02)
#define VSZERO_PLUS   ( ((VSFLOAT) 1E-8) )
#define VSZERO_MINUS  (-VSZERO_PLUS)


#define VSMAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define VSMIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define VSABS(a) ( ( (a) >= 0 ? (a) : -(a) ) )

/* infinite float */
#ifndef inf
#ifdef INFINITY
/* ANSI C99 (gcc -std=c99) */
#define inf INFINITY
#else
/* ANSI C90 (gcc -ansi) */
#define inf DBL_MAX
#endif
#endif

/* -------------------------------------------------------------------------- */
/* VSPstartup version information */
/* -------------------------------------------------------------------------- */
#define VS_DATE "January 1, 2026"
#define VS_MAIN_VERSION 1
#define VS_SUB_VERSION 0
#define VS_SUBSUB_VERSION 0

/* -------------------------------------------------------------------------- */
/* napheap version information */
/* -------------------------------------------------------------------------- */
#define NAPHEAP_DATE "June 10, 2023"
#define NAPHEAP_MAIN_VERSION 5
#define NAPHEAP_SUB_VERSION 0
#define NAPHEAP_SUBSUB_VERSION 0

/* error codes from vs */
#define VS_SUCCESS                                       (0)
#define VS_QR_METHOD_DID_NOT_CONVERGE                    (1)
#define VS_OUT_OF_MEMORY                                 (2)
#define VS_QR_DIAGONALIZATION_FAILS                      (3)
#define VS_NOT_ENOUGH_SPACE_FOR_QR_DIAGONALIZATION       (4)

/* error codes from napheap */
#define NAPHEAP_STATUS_OK                                (0)
#define NAPHEAP_STATUS_INFEASIBLE                      (300)
#define NAPHEAP_STATUS_UNBOUNDED                       (301)
#define NAPHEAP_STATUS_OUT_OF_MEMORY                   (302)
#define NAPHEAP_STATUS_INVALID_N                       (303)
#define NAPHEAP_STATUS_INVALID_BOUNDS                  (304)
#define NAPHEAP_STATUS_INVALID_B                       (305)
#define NAPHEAP_STATUS_INVALID_D                       (306)
#define NAPHEAP_STATUS_CONFLICTING_PARM_D              (307)
#define NAPHEAP_STATUS_MISSING_PARM                    (308)
#define NAPHEAP_STATUS_MISSING_STAT                    (309)

/* other napheap definitions */
#define NAPFLOAT VSFLOAT
#define NAPINT VSINT
#define NONE    0
#define SECANT1 1
#define SECANT2 2
#define VARFIX  3
#define NEWTON  4

#define PRIVATE static
#define NAPHEAP_STATUS_KEEP_LOOKING     (-1)

#define LOOP_OVER_UNFIXED \
nk = nuf ; \
nuf = 0 ; \
for (k = 0; k < nk; k++) \
{ \
    const NAPINT j = uf [k] ;

#define END_UNFIXED_LOOP \
    uf [nuf] = j ; \
    nuf++ ; \
}

#define LOOP_OVER_UNFIXED_NO_UF_UPDATE \
for (k = 0; k < nuf; k++) \
{ \
    const NAPINT j = uf [k] ;


#define END_LOOP \
}

#define FIX_LO \
if ( x [j] == VSZERO ) \
{ \
    continue ; \
}

#define FIX_HI \
if ( x [j] == VSONE ) \
{ \
    b -= aj ; \
    continue ; \
}

#define VAR_FIX_LO \
if ( x [j] == VSZERO ) \
{ \
    const NAPFLOAT t = aj*aj ; \
    den -= t ; \
    num -= t*ay [j] ; \
    continue ; \
}

#define VAR_FIX_HI \
if ( x [j] == VSONE ) \
{ \
    b -= aj ; \
    const NAPFLOAT t = aj*aj ; \
    den -= t ; \
    num -= t*ay [j] ; \
    continue ; \
}

#define SLOPE_UPDATE_LO_NO_DEN \
const NAPFLOAT t = y [j] - lambda*aj ; \
if ( t < VSZERO ) \
{ \
    x [j] = VSZERO ; \
} \
else \
{ \
    if ( t > VSONE ) \
    { \
        x [j] = VSONE ; \
        slope +=aj ; \
    } \
    else \
    { \
        x [j] = t ; \
        slope += aj*t ;

#define SLOPE_UPDATE_HI_NO_DEN \
const NAPFLOAT t = y [j] - lambda*aj ; \
if ( t > VSONE) \
{ \
    x [j] = VSONE ; \
    slope +=aj ; \
} \
else \
{ \
    if ( t < VSZERO ) \
    { \
        x [j] = VSZERO ; \
    } \
    else \
    { \
        x [j] = t ; \
        slope += aj*t ;

#define SLOPE_UPDATE_LO \
SLOPE_UPDATE_LO_NO_DEN \
    } \
}

#define SLOPE_UPDATE_LO_DEN \
SLOPE_UPDATE_LO_NO_DEN \
        den += aj*aj ; \
    } \
}

#define SLOPE_UPDATE_HI \
SLOPE_UPDATE_HI_NO_DEN \
    } \
}

#define SLOPE_UPDATE_HI_DEN \
SLOPE_UPDATE_HI_NO_DEN \
        den += aj*aj ; \
    } \
}

#define BRACKET \
lambdaL = -NAPINF ; \
lambdaR = NAPINF ; \
if ( slope0 < blo ) /* optimal multiplier < 0 */ \
{ \
     b = blo ; \
     slope0 -= b ; \
     dL = NAPINF ; \
     dR = slope0 ; \
     lambdaR = NAPZERO ; \
} \
else if ( slope0 > bhi ) /* optimal multiplier > 0 */ \
{ \
    b = bhi ; \
    slope0 -= b ; \
    dL = slope0 ; \
    dR = -NAPINF ; \
    lambdaL = NAPZERO ; \
} \
else /* lambda = 0 is optimal, blo <= a'x <= bhi */ \
{ \
    lambda = NAPZERO ; \
    for (k = 0; k < n; k++) \
    { \
        const NAPFLOAT t = y [k] ; \
        if ( t < VSZERO ) \
        { \
            x [k] = VSZERO ; \
        } \
        else if ( t > VSONE ) \
        { \
            x [k] = VSONE ; \
        } \
        else \
        { \
            x [k] = t ; \
        } \
    } \
    status = NAPHEAP_STATUS_OK ; \
    return (nap_wrapup (status, napdata, lambda)) ; \
}

/* infinite float */
#define NAPINF inf
#define NAPMAX VSMAX
#define NAPMIN VSMIN
// #define FALSE VSFALSE
// #define TRUE VSTRUE
#define napheap_malloc vs_malloc
#define napheap_free vs_free
#define nap_dot vs_dot
#define NAPZERO VSZERO

/* ========================================================================== */
/* === vs structures ======================================================== */
/* ========================================================================== */

typedef struct VSgivens    /* Givens rotations in diagonalization */
{
    VSFLOAT         *gx ; /* first Givens factor, size=max_its*nalloc */
    VSFLOAT         *gy ; /* second Givens factors, same size as gx */
    short           *gz ; /* 1 or 2, type of Givens transform, size of gx */
    VSFLOAT         *gs ; /* fast Givens scale factors, same size as gx */
    VSINT           *gi ; /* row indices associated with scale factors */
    VSINT           *gj ; /* indexed by QR pass, points into Givens factors*/
    VSINT           *gk ; /* index of last diagonal element in QR iteration */
    VSINT           *gf ; /* index of first diagonal element in QR iteration */
} VSgivens ;

typedef struct VScom
{
    VSparm        *Parm ; /* parameters */
    VSgivens         *G ; /* Givens rotations in diagonalization */
    VSINT        nscale ; /* number of scalings in diagonalization */
    VSINT         npass ; /* number of qr iterations */
    VSINT       ngivens ; /* number of Givens rotations */
    VSINT           *wi ; /* integer work array of size n */
    VSFLOAT         *wx ; /* float work array of size n */
} VScom ;

/* napheap structures */

typedef struct NAPparm_struct
{
    /* T => print status of run */
    int PrintStatus ;

    /* T => print statistics of the run */
    int PrintStat ;

    /* T => print the parameters */
    int PrintParm ;

    /* T => use arrays from a prior run */
    int use_prior_data ;

    /* T => do not free arrays from current run, save them in napheap
            input structure napdata; this is useful when solving the
            same problem but with a different linear cost vector y */
    int return_data ;

    /* F => treat input argument lo as -infinity */
    int loExists ;

    /* F => treat input argument hi as +infinity */
    int hiExists ;

    /* F => no linear constraints bl <= a'x <= bu */
    int Aexists ;

    /* T => diagonal of objective Hessian is all positive, and the
            MATLAB expression all(d>0) is assumed to be true.
       F => there could be zeros on the diagonal of the Hessian */
    int d_is_pos ;

    /* T => diagonal of objective Hessian is zero.
            Note: The d argument of napheap is ignored. It could be NULL
       F => there could be nonzeros on the diagonal of the Hessian */
    int d_is_zero ;

    /* T => diagonal of objective Hessian is identically 1.
            Note: The d argument of napheap is ignored. It could be NULL */
    int d_is_one ;

    /* Upper bound on the number of Newton or variable fixing iterations.
       After performing these K iterations, the code switches to the break
       point search. The Newton algorithm also switches to the breakpoint
       search if two iterates are generated on opposite sides of the root. */
    int K ;

    /* T => use the Newton to update lambda, switch to the break point
       search either at iteration K or at the first time that two
       Newton iterates lie on opposite side of the root (whichever
       happens first)
       F => use K iterations of the variable fixing algorithm before
       switching to a breakpoint search */
    int newton ;

    /* An ordinary Newton iteration often converges monotonically to the
       root. Hence, the Newton iterates are scaled by a factor greater than 1
       in an effort to produce an iterate on the opposite side of the root.
       When such an iterate is generated, the algorithm switches to a
       breakpoint search. The scaled iteration is given by
                    
       lambda_new = lambda_old - newton_scale*L'(lambda_old)/L''(lambda_old) */

    NAPFLOAT newton_scale ;

    /* To combat the effect of rounding errors, recompute a2sum whenever
       a2sum <= decay*a2max, the maximum a2sum encountered */
    double decay ;

     /* T => perform an iterative refinement step to improve solution accurary
        F => no refinement */
    int refine ;

    /* The solution is refined if |a'x - b| > err (only if refine is true) */
    NAPFLOAT err ;

    /* Check the input parameters for dj < 0 or hij < loj or bhi < blo */
    int check ;

} NAPparm ;

/* -------------------------------------------------------------------------- */
/* statistics returned to the user */
/* -------------------------------------------------------------------------- */

typedef struct NAPstat_struct
{
    NAPINT     nkf ;   /* number of known free variables */
    NAPINT   nfree ;   /* number of free variables in initial heap */
    NAPINT  nbound ;   /* number of bound variables in initial heap */
    NAPINT   nbrks ;   /* number break points to reach initial solution */
    NAPINT nrefine ;   /* number of break points during refinement */
    NAPINT nvarfix ;   /* number of variable fixing iterations */
    NAPINT nnewton ;   /* number of Newton iterations performed */
    NAPINT nsecant ;   /* number of secant steps performed */

    /* for error reporting */
    int     status ;   /* returned status from napheap */
    NAPINT  kerror ;   /* invalid n, or index where error found */
    NAPFLOAT lobad ;   /* bad value of lo, if any */
    NAPFLOAT hibad ;   /* bad value of hi, if any */
    NAPFLOAT  dbad ;   /* bad value of d, if any */

} NAPstat ;

/* -------------------------------------------------------------------------- */
/* data structure for input and output to napheap */
/* -------------------------------------------------------------------------- */

typedef struct NAPdata_struc
{   
    NAPFLOAT       *x  ; /* size n. solution (output) */
    NAPFLOAT    lambda ; /* input: starting multiplier guess if parameter
                            use_lambda = TRUE, otherwise ignored on input.
                            output: final multiplier (scalar) */
    NAPINT           n ; /* problem dimension */
    NAPFLOAT        *c ; /* size n. linear term in objective function */
    NAPFLOAT        *d ; /* size n. diagonal of objective Hessian */
    NAPFLOAT        *a ; /* size n. linear constraint vector */
    NAPFLOAT       blo ; /* lower bound for a'x */
    NAPFLOAT       bhi ; /* upper bound for a'x */
    NAPFLOAT       *lo ; /* size n. lower bounds for x */
    NAPFLOAT       *hi ; /* size n. upper bounds for x */
    NAPstat      *Stat ; /* structure for statistics */
    NAPparm      *Parm ; /* structure for parameters */
    NAPFLOAT    *xWork ; /* NULL => let code allocate NAPFLOAT memory
                            otherwise use xWork for NAPFLOAT memory.
                            The amount of memory needed depends on the
                            parameters in the Parm structure. See the
                            memory allocation section of the code for the
                            precise memory requirement (max size 5n). */
    NAPINT      *iWork ; /* NULL => let code allocate NAPINT memory
                            otherwise use iWork for NAPINT memory (4*n+1).*/

    /* the following are used internally when the code allocates memory */
    NAPFLOAT     *x_created ;
    NAPFLOAT     *c_created ;
    NAPFLOAT *xWork_created ;
    NAPINT   *iWork_created ;
    /* sum ak*ak/dk, k = 1:n, from a previous run */
    NAPFLOAT         akakdk ;
} NAPdata ;

/* Ovi's structures */

typedef struct
{
    VSINT *queue ;
    VSINT  size_max ;  // maximum number of members in the queue
    VSINT  ptr_start ; // index in queue of oldest member. Initialized as 0
    VSINT  ptr_end ;   // index in queue of newest member. Initialized as -1 
} Queue ;

typedef struct
{
    VSINT   m ;    // number of rows
    VSINT   n ;    // number of columns
    VSLINT *p ;    // column pointers (size n+1)
    VSINT  *i ;    // row indices, size nz
    VSLINT nz ;    // number of nonzeros
} CSC ;


void vs_print_parms
(
    VSparm *Parm /* VSparm structure to be printed */
) ;

void * vs_malloc
(
    int *status,
    VSINT     n,
    int    size
) ;

void * vs_lmalloc
(
    int *status,
    VSLINT    n,
    int    size
) ;

void * vs_Lmalloc
(
    int *status,
    int64_t   n,
    int    size
) ;

void vs_free
(
    void * p
) ;

void vs_printAMATLAB
(
    VSINT   const ncol, /* number of cols in A */
    VSLINT  const  *Ap, /* size ncol+1, column pointers */
    VSINT   const  *Ai, /* size Ap [ncol], row indices for A */
    VSFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
) ;

VSFLOAT vs_max
(
    VSFLOAT const *x, /* VSFLOAT vector */
    VSINT   const  n  /* VSINT length of vector */
) ;

VSFLOAT vs_norm
(
    VSFLOAT const *x, /* return Euclidean norm of x */
    VSINT   const  n  /* dimension */
) ;

void vs_initx
(
    VSFLOAT      *x,  /* array to be initialized */
    VSFLOAT const s,  /* scalar */
    VSINT   const n   /* length of x */
) ;

void vs_initi
(
    VSINT      *x,  /* array to be initialized */
    VSINT const s,  /* scalar */
    VSINT const n   /* length of x */
) ;

void vs_initLx
(
    VSFLOAT      *x,  /* array to be initialized */
    VSFLOAT const s,  /* scalar */
    VSLINT  const n   /* length of x */
) ;

void vs_copyx
(
    VSFLOAT       *x, /* output of copy */
    VSFLOAT const *y, /* input of copy */
    VSINT   const  n  /* length of vectors */
) ;

void vs_vecadd
(
     VSFLOAT       *x,
     VSFLOAT       *y,
     VSFLOAT const *b,
     VSINT   const  n
) ;

void vs_scale
(
    VSFLOAT       *x,  /* scaled array */
    VSFLOAT const  s,  /* scale */
    VSINT   const  n   /* length of x */
) ;

void vs_project
(
    VSFLOAT       *x, /* scaled array */
    VSFLOAT const *a, /* project onto orthogonal complement of a */
    VSFLOAT const aa, /* 1/a'a */
    int     const a_is_one,
    VSINT   const  n  /* dimension */
) ;

VSFLOAT vs_dot
(
    VSFLOAT const *a,
    VSFLOAT const *x,
    VSINT   const  n
) ;

void vs_translate
(   
    VSFLOAT    *new,
    VSFLOAT const t,
    VSINT   const n
) ;

void vs_boundary /* return original cost and its improvement */
(
    VSFLOAT           *x, /* feasible point for (GB), output described above */
    VSINT            *Ai, /* row numbers */
    VSLINT           *Ap, /* column pointers */
    VSFLOAT          *Ax, /* edge weights,   NULL => all entries in Ax are 1 */
    VSFLOAT           *a, /* vertex weights, NULL => all entries in a are 1 */
    VSINT              n, /* problem dimension */
    VSFLOAT         *A12, /* linear term in objective of (GB) */
    VSFLOAT           *g, /* work array for storing gradient vector */
    int const PrintLevel
) ;

VSINT vs_interior
(
   VSINT    ilower, /* search x for indices > ilower and < n */
   VSFLOAT      *x, /* search for 0 < x [i] < 1 */
   VSINT         n  /* dimension of x */
) ;

void vs_house
(
    VSFLOAT *A, /* input:  n-by-n dense matrix
                    output: tridiagonalization plus Household product in A */
    VSFLOAT *d, /* n-by-1 vector (output) */
    VSFLOAT *u, /* n-by-1 vector (output) */
    VSINT    n, /* matrix dimension */
    VScom  *com  /* common data */
) ;

int vs_diag /* return: 0 (convergence tolerance satisfied)
                              -5 (algorithm did not converge)
                              -6 (insufficient space allocated for Givens) */
(
    VSFLOAT   *d, /* diagonal of tridiagonal matrix, output: eigenvalues */
    VSFLOAT   *u, /* superdiagonal of tridiagonal matrix, use u (0:n-2) */
    VSINT      n, /* dimensional of tridiagonal matrix */
    VScom   *com  /* common data */
) ;

void vs_DenseMult
(
    VSFLOAT *y,     /* m by 1 product, output */
    VSFLOAT *x,     /* n by 1 given vector */  
    VSFLOAT *V,     /* dense m by n matrix */
    VSINT    m,     /* number of rows */
    VSINT    n      /* number of columns */
) ;

void vs_tDenseMult
(
    VSFLOAT *y,     /* n by 1 product, output */
    VSFLOAT *x,     /* m by 1 given vector */
    VSFLOAT *V,     /* dense m by n matrix */
    VSINT    m,     /* number of rows */
    VSINT    n      /* number of columns */
) ;

void vs_GivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT       n  /* dimension of matrix */
) ;

void vs_tGivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT  nscale, /* number of scaling operations performed */
    VSINT       n  /* dimension of matrix */
) ;

void vs_Amult
(
    VSFLOAT     *Av, // Stores the result, A * v
    VSFLOAT      *v, // Vector being multiplied by A
    VSLINT      *Ap, // Row pointers
    VSINT       *Ai, // Column indices
    VSFLOAT     *Ax, // Matrix entries 
    VSINT         n  // Number of columns
) ;

void vs_A1mult
(       
    VSFLOAT     *A1, // Stores the product
    VSLINT      *Ap, // Row pointers
    VSFLOAT     *Ax, // Matrix entries
    VSINT         n, // Number of columns
    int    E_is_one  // Flag (set to 1 for binary matrix, 0 otherwise) 
) ;

int vs_sort_cols
(
    VSLINT      *Ap, /* given column pointers for symmetric matrix */
    VSINT       *Ai, /* sorted row indices (spaced malloc'd by user) */
    VSINT       *Ax, /* sorted numerical entries (spaced malloc'd by user) */
    VSINT const *Bi, /* given row indices (unsorted) */
    VSINT const *Bx, /* given entries (unsorted), Bx = NULL => Ax = NULL */
    VSINT         n  /* number of rows in matrix */
) ;

void vs_daxpy
(
    VSFLOAT       *y,
    VSFLOAT const *x,
    VSFLOAT const  a,
    VSINT   const  n
)  ;

void vs_GivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT       n  /* dimension of matrix */
) ;

void vs_tGivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT  nscale, /* number of scaling operations performed */
    VSINT       n  /* dimension of matrix */
) ;

void vs_printA
(
    VSINT   const ncol, /* number of cols in A */
    VSLINT  const  *Ap, /* size ncol+1, column pointers */
    VSINT   const  *Ai, /* size Ap [ncol], row indices for A */
    VSFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    int     const diag, /* 1 = add diag */
    char    const *what /* name of the matrix */
) ;

void vs_printxMATLAB
(
    VSFLOAT *x, /* numerical entries in the vector */
    VSINT    n, /* dimension of the vector */
    char *what  /* name of the vector */
) ;

#define TIME( A ) clock_gettime(CLOCK_MONOTONIC, &(A)) ;

// Input:
// Ai, Ap, n - the 
// part[i] is the partition (0 or 1 2) vertex i is in
// Ai are the row indices of the adjacency matrix
// Ap holds the column pointers of the adjacency matrix
// n is the number of vertices of the graph 
void vs_separator
(
    int                *part, // array with entries 0, 1, or 2
    VSINT                *Ai, // Original graph
    VSLINT               *Ap, // Original graph
    VSINT                  n, // Original graph n
    VSINT          *vtx_wgts,
    int      try_both_covers, // If nonzero: compute both Konig covers (BFS
                              // from unmatched cols and from unmatched rows)
                              // and pick the one with smallest weight,
                              // tiebreaking on smallest imbalance.
                              // If zero: original logic, single cover from
                              // BFS on unmatched columns.
    int            get_stats,  // Binary flag to determine if (somewhat
                               // expensive) stats should be computed

    // Ignored if get_stats is false. In particular, they can be NULL
    VSINT             *n_bip,
    VSLINT           *nz_bip,
    VSINT         *wgt_cover,
    VSINT             *imbal,
    VSFLOAT *time_entire_ret
) ;

void vs_separator_approx
(
    int          *part, // array with entries 0, 1, or 2
    VSINT          *Ai, // Original graph
    VSLINT         *Ap, // Original graph
    VSINT            n, // Original graph n
    VSINT    *vtx_wgts,
    int       get_stats,  // Binary flag to determine if (somewhat expensive)
                          // stats should be computed

    // Ignored if get_stats is false. In particular, they can be NULL
    VSINT              *num_vert_bip,
    VSLINT                   *nz_bip,
    VSINT                 *wgt_cover,
    VSINT                     *imbal,
    VSFLOAT         *time_entire_ret
) ;

void vs_print_TF
(
    int TF /* TRUE or FALSE (or EMPTY if not True or False) */
) ;

/* prototypes for napheap */

int napheap  /* napheap_print_status (napdata) gives the status of the run */
(
    NAPdata *napdata
) ;

NAPdata * napheap_setup (void) ;

void napheap_terminate
(
    NAPdata **DataHandle
) ;

/* set default parameter values */
void napheap_default
(
    NAPparm *Parm
) ;

void napheap_print_status
(
    NAPdata *Data
) ;

/* print the napheap statistics */
void napheap_print_stat
(
    NAPdata *Data
) ;

/* print the parameters */
void napheap_print_parm
(
    NAPdata *Data
) ;

PRIVATE int nap_wrapup
(
    int         status,
    NAPdata   *napdata,
    NAPFLOAT    lambda
) ;

PRIVATE int napsearch   /* return status of solution process */
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT   *lambda, /* optimal multiplier (input = guess, output = sol)*/
    NAPFLOAT   lambdaL, /* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT   lambdaR, /* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT         B, /* a'x = B */
    NAPFLOAT     slope, /* derivative of dual function */
    NAPFLOAT const *ay, /* size n. linear term in objective function yj/aj */
    NAPFLOAT const *ad, /* size n, diagonal of objective Hessian aj/dj */
    NAPFLOAT const  *a, /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPINT   const   n, /* problem dimension */
    NAPFLOAT *breakpts, /* size n, breakpts [j] is br_hi[j] or br_lo[j] */
    NAPFLOAT *breaknext,/* size n, TRUE if 2 break points in search dir. */
    NAPINT         nuf, /* number of unfixed variables */
    NAPINT *known_free, /* location of first known free variable in uf */
    NAPINT     *n_free, /* number of indices in free heap */
    NAPINT    *n_bound, /* number of indices in bound heap */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT         *kf, /* size n, list of known free variables */
    NAPINT  *free_heap, /* size n+1, heap of free variables */
    NAPINT *bound_heap, /* size n+1, heap of bound variables */
    NAPstat      *Stat, /* solution statistics */
    NAPparm      *Parm, /* parameters */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists  /* TRUE when there are upper bound for x */
) ;

PRIVATE void napsearch0
(
    NAPFLOAT              *x, /* size n. solution (output) */
    NAPFLOAT         *lambda, /* optimal multiplier (in = guess, out = sol) */
    NAPFLOAT         lambdaL, /* optimal multiplier lies on [lambdaL, lambdaR]*/
    NAPFLOAT         lambdaR, /* optimal multiplier lies on [lambdaL, lambdaR]*/
    NAPFLOAT              *b, /* a'x = b */
    NAPFLOAT          *slope, /* slope (exclude break point term) */
    NAPFLOAT        *slopelo, /* lower subdifferential associated with break */
    NAPFLOAT        *slopehi, /* upper subdifferential associated with break */
    NAPINT   const         n, /* problem dimension */
    NAPFLOAT const        *a, /* size n, linear constraint vector */
    NAPFLOAT const       *lo, /* size n, lower bounds on x */
    NAPFLOAT const       *hi, /* size n, upper bounds on x */
    NAPFLOAT const *breakpts, /* size n, breakpts [j] is br_hi[j] or br_lo[j] */
    NAPINT               *uf, /* size n, list of unfixed variables */
    NAPINT         *nunfixed, /* number of unfixed variables */
    NAPINT       *bound_heap, /* size n+1, heap of bound variables */
    int              goright, /* if true, move right, else move left */
    NAPstat            *Stat, /* solution statistics */
    int      recompute_slope, /* if true, then recompute slope */
    NAPparm            *Parm, /* parameters */
    int      const  loExists, /* TRUE when there are lower bound for x */
    int      const  hiExists  /* TRUE when there are upper bound for x */

) ;

PRIVATE void napsearch1
(
    NAPFLOAT           *x,/* size n. solution (output) */
    NAPFLOAT      *lambda,/* optimal multiplier (input = guess, output = sol) */
    NAPFLOAT      lambdaL,/* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT      lambdaR,/* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT           *b,/* a'x = b */
    NAPFLOAT       *slope,/*  slope excluding term associated with break point*/
    NAPFLOAT     *slopelo,/* lower subdifferential associated with break point*/
    NAPFLOAT     *slopehi,/* upper subdifferential associated with break point*/
    NAPFLOAT const     *y,/* size n. linear term in objective function */
    NAPFLOAT const     *d,/* size n, diagonal of objective Hessian */
    NAPFLOAT const     *a,/* size n, linear constraint vector */
    NAPFLOAT const    *lo,/* size n, lower bounds on x */
    NAPFLOAT const    *hi,/* size n, upper bounds on x */
    NAPINT   const      n,/* problem dimension */
    NAPINT              K,/* number iterations in variable fixing method */
    NAPFLOAT    *breakpts,/* size n, breakpts [j] is br_hi[j] or br_lo[j] */
    NAPINT            *uf,/* size n, list of unfixed variables */
    NAPINT      *nunfixed,/* number of unfixed variables */
    NAPINT     *free_heap,/* size n+1, heap of free variables */
    NAPINT    *bound_heap,/* size n+1, heap of bound variables */
    int           goright,/* if true, move right, else move left */
    NAPstat         *Stat,/* solution statistics */
    NAPFLOAT const *br_lo,/* size n, br_lo[j] lo break pt for x[j] */
    NAPFLOAT const *br_hi,/* size n. br_hi[j] hi break pt for x[j] */
    int   recompute_slope,/* if true, then recompute slope */
    NAPparm         *Parm,/* parameters */
    int   const  loExists,/* TRUE when there are lower bound for x */
    int   const  hiExists /* TRUE when there are upper bound for x */
) ;

PRIVATE void napfix0
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT const  *a, /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT        *b, /* b - ai*xi for fixed variables */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT   *nunfixed, /* number of unfixed variables */
    NAPFLOAT const *breakpts /* break points */
) ;

PRIVATE void napfix1
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT const  *a, /* size n, linear constraint vector */
    NAPFLOAT const  *y, /* size n. linear term in objective function */
    NAPFLOAT const  *d, /* size n, diagonal of objective Hessian */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT        *b, /* b - ai*xi for fixed variables */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT   *nunfixed, /* number of unfixed variables */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists, /* TRUE when there are upper bound for x */
    NAPFLOAT const *breakpts /* break points */
) ;

PRIVATE void napsolution
(
    NAPFLOAT        *x, /* size n, solution (output) */
    NAPFLOAT const  *a, /* size n, coefficient of linear constraint */
    NAPFLOAT const  *y, /* size n, linear term in objective */
    NAPFLOAT const *ay, /* size n, linear term y in objective over a */
    NAPFLOAT const *ad, /* size n, coefficient of linear constraint over d */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPINT   const *uf, /* size n, unfixed variables */
    NAPINT   const *kf, /* size n, known free variables */
    NAPINT         nuf, /* number of unfixed variables */
    NAPINT         nkf, /* number of known free variables */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists, /* TRUE when there are upper bound for x */
    int const d_is_one  /* TRUE when diagonal is identically one */
) ;

PRIVATE void napsolution0
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT const  *a, /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT         B, /* the constraint is A'X = B */
    NAPFLOAT *breakpts, /* size n, breakpts [j] is br_hi[j] or br_lo[j] */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT         nuf, /* number of unfixed variables */
    NAPINT  *free_heap, /* size n+1, heap of free variables */
    NAPINT *bound_heap, /* size n+1, heap of bound variables */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists  /* TRUE when there are upper bound for x */
) ;

PRIVATE void napsolution1
(
    NAPFLOAT        *x,  /* size n. solution (output) */
    NAPFLOAT const  *y,  /* size n. linear term in objective function */
    NAPFLOAT const  *d,  /* size n, diagonal of objective Hessian */
    NAPFLOAT const  *a,  /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPFLOAT     slope, /* slope of dual function determines search direction */
    NAPINT   const *kf,/* size n, known free variables */
    NAPINT         nkf, /* number of known free variables */
    NAPINT  *free_heap,
    NAPINT *bound_heap,
    NAPINT      n_free,
    NAPINT     n_bound,
    int const d_is_one
) ;

PRIVATE void napminheap_build
(
    NAPINT      *heap,  /* on input, unsorted set of indices.  size nheap+1 */
    NAPFLOAT const *X,  /* numbers to sort */
    NAPINT      nheap   /* number of elements to build into the heap */
) ;

PRIVATE NAPINT napminheap_delete  /* return new size of heap */
(
    NAPINT      *heap,  /* indices into X, 1..n on input.  size nheap+1 */
    NAPFLOAT const *X,  /* not modified */
    NAPINT      nheap   /* number of items in heap */
) ;

PRIVATE NAPINT napminheap_add
(
    NAPINT      *heap,  /* size n, containing indices into X.  size nheap+2 */
    NAPFLOAT const *X,  /* not modified */
    NAPINT       leaf,  /* the new leaf */
    NAPINT      nheap   /* number of elements in heap not counting new one */
) ;

PRIVATE void napminheapify
(
    NAPINT      *heap,  /* size nheap+1, containing indices into X */
    NAPFLOAT const *X,  /* not modified */
    NAPINT          p,  /* start at node p in the heap */
    NAPINT      nheap   /* heap [1 ... n] is in use */
) ;

PRIVATE void napmaxheap_build
(
    NAPINT      *heap,  /* on input, unsorted set of elements. size nheap+1 */
    NAPFLOAT const *X,  /* not modified */
    NAPINT      nheap   /* number of elements to build into the heap */
) ;

PRIVATE NAPINT napmaxheap_delete  /* return new size of heap */
(
    NAPINT      *heap,  /* containing indices into X, size nheap+1 */
    NAPFLOAT const *X,  /* not modified */
    NAPINT      nheap   /* number of items in heap */
) ;

PRIVATE NAPINT napmaxheap_add
(
    NAPINT      *heap,  /* size nheap+2, containing indices into X */
    NAPFLOAT const *X,  /* not modified */
    NAPINT       leaf,  /* the new leaf */
    NAPINT      nheap   /* number of elements in heap not counting new one */
) ;

PRIVATE void napmaxheapify
(
    NAPINT      *heap,  /* size nheap+1, containing indices into X */
    NAPFLOAT const *X,  /* not modified */
    NAPINT          p,  /* start at node p in the heap */
    NAPINT      nheap   /* heap [1 ... nheap] is in use */
) ;

PRIVATE void nap_infeasible
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT const  *a, /* size n. linear constraint vector */
    NAPFLOAT const *lo, /* size n. lower bounds for x */
    NAPFLOAT const *hi, /* size n. upper bounds for x */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT         nuf, /* number of unfixed variables */
    NAPFLOAT     minax,
    NAPFLOAT     maxax,
    NAPFLOAT       blo, /* lower bound for a'x */
    NAPFLOAT       bhi  /* upper bound for a'x */
) ;

PRIVATE int nap_work
(
    NAPINT  const      n,
    NAPFLOAT      *xWork,     /* user-provide workspace */
    NAPFLOAT       **xw2,     /* return malloc'ed space */

    /* the 6 work arrays, each of size n */
    NAPFLOAT        **ad,
    NAPFLOAT        **ay,
    NAPFLOAT     **br_hi,
    NAPFLOAT     **br_lo,
    NAPFLOAT  **breakpts,
    NAPFLOAT **breaknext
) ;

NAPFLOAT napheap_check /* returns largest of the 3 error values */
(   
    NAPFLOAT     *errb, /* relative error in b (only return if not NULL) */
    NAPFLOAT     *errB, /* relative error in bounds (only return if not NULL) */
    NAPFLOAT     *errc, /* relative error in c (only return if not NULL) */
    NAPdata   *napdata
) ;

/* prototypes for codes to evaluate a minimum vertex cover */

static VSINT *Randperm (VSINT n) ;

static void Accumulate (VSLINT *p, VSLINT *c, VSINT n) ;

static void Elapsed_Time
(
    struct timespec start, struct timespec end, VSFLOAT *elapsed
) ;

static Queue *Queue_Init(VSINT size_max) ;

static void Queue_Destroy(Queue *Queue) ;

static void Queue_Push(Queue *Queue, VSINT memb) ;

static VSINT Queue_Pop(Queue *Queue) ;

static CSC *CSC_Init(VSINT m, VSINT n, VSLINT nz) ;

static void CSC_Destroy(CSC *A) ;

static CSC *CSC_Transpose (const CSC *A) ;

static CSC *CSC_Extract_Bipartite_Graph
(
    VSINT              *part,
    CSC                    *A,
    VSINT **localBipLabelRet
) ;

void MinHeap_Bubble_Down_v0
(
    VSINT  *heap, 
    VSINT  *lookup,
    VSINT   numNodes,
    VSINT  *c,
    VSINT   p
) ;

void MinHeap_Heapify_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT   numNodes,
    VSINT  *c
) ;

void MinHeap_Insert_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT   toInsert
) ;

void MinHeap_Delete_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT   toDel
) ;

void MinHeap_See_Top_Val_v0
(
    VSINT  *heap,
    VSINT  *c,
    VSINT   numNodes,
    VSINT  *topVal
) ;

void MinHeap_Get_Top_Val_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT  *topVal
) ;

void MinHeap_See_Top_Key_v0
(
    VSINT  *heap,
    VSINT  *c,
    VSINT   numNodes,
    VSINT  *topKey
) ;

void MinHeap_Get_Top_Key_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT  *topKey
) ;

static void Augment (VSINT k, const CSC *A, VSINT *jmatch, VSLINT *cheap,
        VSINT *w, VSINT *js, VSINT *is, VSLINT *ps) ;

static VSINT *Maximum_Matching_Get(const CSC *A, VSINT seed) ;

static int Are_All_Neighbors_In_Cover
(
    CSC             *B,
    VSINT     vtx_col,
    VSINT      *cover  // cover[i] = 1 if row vertex i is in the cover  
) ;

static void Check_Bip_Vertex_Cover(CSC *B, VSINT *cover) ;

static void Check_Matching_Size_Against_Cover_Size
(   
    CSC *B,
    VSINT *matching,
    VSINT *cover
) ;

static void Separator_Check
(
    VSINT *part, // entries are 0, 1, 2. Part 2 is candidate separator
    CSC       *A  // Adjacency matrix
) ;

static void Check_Its_A_Cover(CSC *B, VSINT *cover) ;

void Pricing_Algorithm_Vertex_Cover_Get
(
    CSC            *B,
    VSINT         *w, 
    VSINT     *cover
) ;


int vs_separator_min_weight
(
    int                *part, // array with entries 0, 1, or 2
    VSINT                *Ai, // Original graph
    VSLINT               *Ap, // Original graph
    VSINT                  n, // Original graph n
    VSINT          *vtx_wgts
) ;

#endif

// // IWYU pragma: keep

// #ifndef VSP_IOH
// #define VSP_IOH


// #include "shared_api.h"
// #include <limits.h>
// #include <stddef.h>
// #ifdef MATLAB_MEX_FILE
// #include "mex.h" 
// #endif

// /* same as Int in Mongoose */

// #define Int long

// /* CSparse Macros */
// #ifndef CS_CSC
// #define CS_CSC(A) ((A) && ((A)->nz == -1))
// #endif
// #ifndef CS_TRIPLET
// #define CS_TRIPLET(A) ((A) && ((A)->nz >= 0))
// #endif

// /* --- primary CSparse routines and data structures ------------------------- */
// typedef struct cs_sparse /* matrix in compressed-column or triplet form */
// {
//     csi nzmax; /* maximum number of entries */
//     csi m;     /* number of rows */
//     csi n;     /* number of columns */
//     csi *p;    /* column pointers (size n+1) or col indices (size nzmax) */
//     csi *i;    /* row indices, size nzmax */
//     double *x; /* numerical values, size nzmax */
//     csi nz;    /* # of entries in triplet matrix, -1 for compressed-col */
// } cs;

// cs *cs_add(const cs *A, const cs *B, double alpha, double beta);
// cs *cs_transpose(const cs *A, csi values);

// cs *cs_compress(const cs *T);

// cs *cs_spalloc(csi m, csi n, csi nzmax, csi values, csi triplet);
// cs *cs_spfree(cs *A);

// cs *sanitizeMatrix(cs *compressed_A, bool symmetricTriangular,
//                    bool makeEdgeWeightsBinary);
// void removeDiagonal(cs *A);
// // Requires A to be a triangular matrix with no diagonal.
// cs *mirrorTriangular(cs *A);

// /* My own functions */
// VSP_problem* VSP_read_MM_file(const char* filename, VSP_options* options);

// #endif
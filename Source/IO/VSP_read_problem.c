#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared_api.h"
#include "io_internal.h"
#include "VSP_read_problem.h"




VSP_problem* VSP_read_problem(char* file_name, VSP_options* options) {

    char* dot = strrchr( file_name, '.' ) ;
    VSP_problem *problem = NULL ;
    if ( dot [1] == 'm' && dot [2] == 't' && dot [3] == 'x' )
    {
        problem = VSP_read_MM_file(file_name, options) ;

        // Clean up and write file name
        char *lastSlash = strrchr(file_name, '/') ;
        strncpy(problem->name, lastSlash + 1, (dot - lastSlash) ) ;
        (problem->name)[dot - lastSlash - 1] = '\0' ;
        return problem ;
    }
    else if( dot [1] == 'o' && dot [2] == 'v' && dot [3] == 'i' )
    {
        problem = Read_Dot_Ovi_File(file_name, options) ;

        // Clean up and write file name
        char *lastSlash = strrchr(file_name, '/') ;
        strncpy(problem->name, lastSlash + 1, (dot - lastSlash) ) ;
        (problem->name)[dot - lastSlash - 1] = '\0' ;
        return problem ;
    }
    else
    {
        printf( "\n\n ERROR: File format not recognized by VSP solver.\n\n" ) ;
        exit(1) ;
    }
}

/*------------------------------------------------------------------------*/
// .ovi files have the following format (copied from can___61):
// n: 61    nz: 496
// Ap:
// 0
// 5
// 13
// 21
// 45
// 58
// ...
// 487
// 492
// 496


// Ai:
// 13
// 15
// 19
// 23
// 24
// ...
// 55
// 56

/*------------------------------------------------------------------------*/

VSP_problem* VSP_read_MM_file(const char* filename, VSP_options* options)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count += 1 ;
    #endif

    MM_typecode matcode;

    /* Open the file */
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Cannot read file\n");
        exit(1) ;
    }

    /* Read MM banner */
    if (mm_read_banner(file, &matcode) != 0) {
        printf("Could not process Matrix Market banner\n");
        exit(1) ;
    }

    /* Check that given matrix is of the correct format */
    if (!mm_is_matrix(matcode))
    {
        printf("Error: File not a matrix:  %s\n", filename);
        exit(1);
    }
    if (!mm_is_sparse(matcode))
    {
        printf("Error: Matrix is not sparse:  %s\n", filename);
        exit(1);
    }
    if (mm_is_complex(matcode))
    {
        printf("Error: Matrix is complex, must be real:  %s\n", filename);
        exit(1);
    }


    /* Read matrix size and nonzeroes */
    long M, N, nz;
    if ((mm_read_mtx_crd_size(file, &M, &N, &nz)) != 0)
    {
        printf("Could not parse matrix dimension and size.\n");
        printf("\nFile name: %s", filename) ;
        exit(1) ;
    }
    if (M != N) {
        printf("\nMatrix must be square.");
        printf("\nFile name: %s", filename) ;
        exit(1) ;
    }


    /* Read matrix entries */
    Int *I = (Int *)SuiteSparse_malloc((size_t) nz, sizeof(Int));
    Int *J = (Int *)SuiteSparse_malloc((size_t) nz, sizeof(Int));
    double *val = (double *)SuiteSparse_malloc((size_t) nz, sizeof(double));

    if (!I || !J || !val)
    {
        printf("Error: Ran out of memory");
        SuiteSparse_free(I);
        SuiteSparse_free(J);
        SuiteSparse_free(val);

    }
    mm_read_mtx_crd_data(file, M, N, nz, (long*)I, (long*)J, val, matcode);

    /* Close the file */
    fclose(file);

    /* MM matrix indices are 1-based. Make them 0-based */
    for (Int k = 0; k < nz; k++)
    {
        --I[k];
        --J[k];
        if (mm_is_pattern(matcode))
            val[k] = 1;
    }

    /* Put the above info in a cs matrix A */
    cs *A = (cs *)SuiteSparse_malloc(1, sizeof(cs));
    if (!A) {
        printf("Ran out of memory.\n");
        SuiteSparse_free(I);
        SuiteSparse_free(J);
        SuiteSparse_free(val);
    }

    A->nzmax = nz;
    A->m     = M;
    A->n     = N;
    A->p     = J;
    A->i     = I;
    A->x     = val;
    A->nz    = nz;

    /* Put the matrix in compressed column format */
    cs *compressed_A = cs_compress(A);
    cs_spfree(A);

    if (!compressed_A) {
        printf("Ran out of memory.\n");
    }


    #ifdef VSP_DEBUG
        VSPINT n_dbg = compressed_A->n ;
        if (compressed_A->p[n_dbg] != compressed_A->nz) {
            VSP_error("Nz doesn't match Ap[n]", __FILE__, __LINE__, count) ;
        }
    #endif


    /* ====================================================================== */
    /* ============ September 19, 2023 ====================================== */
    /* ====================================================================== */

    /* This block of code is for dealing with non-symmetric matrices A, which
       arise from directed graphs. What we do is turn the adjacency of the 
       directed graph into the matrix of the corresponding undirected graph
       as follows:
       1) Get the transpose of A
       2) Do A + A^T
       3) Ignore the entries in A + A^T, and just keep the pattern of 1's
    */

    if ( matcode[3] == 'G' ) {
        cs* compressed_A_transpose = cs_transpose(compressed_A, TRUE) ;
        cs* compressed_A_symmetrized = cs_add(
                                            compressed_A,
                                            compressed_A_transpose,
                                            1,
                                            1) ;
        
        VSPINT num_vertices = compressed_A_symmetrized->n ;
        compressed_A_symmetrized->nz = 
                                    compressed_A_symmetrized->p[num_vertices] ;

        removeDiagonal(compressed_A_symmetrized) ;

        VSP_problem* problem =  VSP_construct_problem
                                (
                                    compressed_A_symmetrized->n,
                                    compressed_A_symmetrized->nz,
                                    NULL, NULL, NULL,
                                    options
                                ) ;
        

        for(VSPINT i = 0; i < problem->Graph->n + 1; i++){
            problem->Graph->Ap[i] = (VSPLONG) compressed_A_symmetrized->p[i];
        }

        for(VSPINT i = 0; i < problem->Graph->nz; i++){
            problem->Graph->Ai[i] = (VSPLONG) compressed_A_symmetrized->i[i];
        }

        /* Free memory */
        cs_spfree(compressed_A) ;
        cs_spfree(compressed_A_transpose) ;
        cs_spfree(compressed_A_symmetrized) ;

        problem->Graph->nz = problem->Graph->Ap [problem->Graph->n] ;
        return problem ;
    }
   



    /* ====================================================================== */
    /* ====================================================================== */
    /* ====================================================================== */


    /* Sanitize the matrix */
    cs *sanitized_A = sanitizeMatrix(compressed_A, mm_is_symmetric(matcode), 0);

    if (!sanitized_A) {
        printf("Ran out of memory.");
    }


    VSP_problem* problem = VSP_construct_problem
                            ((VSPINT) sanitized_A->n, (VSPINT) sanitized_A->nz,
                             NULL, NULL, NULL, options);

    for(VSPINT i = 0; i < problem->Graph->n + 1; i++){
        problem->Graph->Ap[i] = (VSPINT) sanitized_A->p[i];
    }

    for(VSPINT i = 0; i < problem->Graph->nz; i++){
        problem->Graph->Ai[i] = (VSPINT) sanitized_A->i[i];
    }

    /* Free memory */
    cs_spfree(compressed_A) ;
    cs_spfree(sanitized_A) ;

    problem->Graph->nz = problem->Graph->Ap [problem->Graph->n] ;
    return problem;

}


// Function to extract base filename
// Parameters:
// - file_path: The full path to the file.
// - dest: Destination buffer to store the base filename.
// - dest_size: Size of the destination buffer.
// Returns:
// - 0 on success, -1 on failure.
VSPINT extract_base_filename(const char *file_path, char *dest, size_t dest_size) {
    if (file_path == NULL || dest == NULL) {
        return -1; // Invalid arguments
    }

    const char *last_slash = strrchr(file_path, '/');
    const char *last_dot = strrchr(file_path, '.');

    // Determine the start of the base filename
    const char *start = (last_slash != NULL) ? last_slash + 1 : file_path;

    // Determine the end of the base filename
    const char *end;
    if (last_dot != NULL && last_dot > start) {
        end = last_dot;
    } else {
        end = file_path + strlen(file_path); // No extension found
    }

    // Calculate the length of the base filename
    size_t len = end - start;

    // Handle edge case where the file starts with a dot (hidden files without extension)
    if (len == 0 && last_dot != NULL && last_dot == start) {
        // The file is something like ".bashrc"
        len = strlen(start);
        end = file_path + strlen(file_path);
    }

    // Ensure we don't exceed the destination buffer
    if (len >= dest_size) {
        len = dest_size - 1; // Reserve space for null terminator
    }

    // Copy the substring into dest
    strncpy(dest, start, len);
    dest[len] = '\0'; // Null-terminate the string

    return 0; // Success
}
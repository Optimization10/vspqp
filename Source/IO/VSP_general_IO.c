#include "shared_api.h"
#include "shared_types.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h>

VSP_solution_stats* VSP_compute_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    VSP_problem *problem,
    VSPFLOAT        time
)
{
    /* Local vars */
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT num_nodes_x = 0 ;
    VSPINT num_nodes_y = 0 ;

    /* Allocate memory */
    VSP_solution_stats *stats = calloc(1, sizeof(*stats)) ;

    /* Compute number of nodes in x and y */
    for(i = 0 ; i < n ; i++) {
        num_nodes_x += x[i] ;
        num_nodes_y += y[i] ;
    }

    /* Fill out number of nodes in stats struct */
    stats->num_nodes_total = n ;
    stats->num_nodes_x     = num_nodes_x ;
    stats->num_nodes_y     = num_nodes_y ;
    stats->num_nodes_sep   = n - num_nodes_x - num_nodes_y ;

    stats->num_edges_total = problem->Graph->nz / 2 ;

    /* Fill out percentages in stats struct */
    stats->percent_nodes_x   = (VSPFLOAT) num_nodes_x * 100 / n ;
    stats->percent_nodes_y   = (VSPFLOAT) num_nodes_y * 100 / n ;
    stats->percent_nodes_sep = (VSPFLOAT) stats->num_nodes_sep * 100 / n ;

    /* Fill out the time */
    stats->time = time ;

    return stats ;
}



void VSP_reset_file(char* name)
{
    FILE* file = fopen(name, "w") ;
    fclose(file) ;
    return ;
}

void VSP_reset_output_files()
{
    VSP_reset_file("output/VSP_results.txt") ;
    VSP_reset_file("output/VSP_error_log.txt") ;
    VSP_reset_file("output/VSP_stats.txt") ;
    return ;
}

/* Dump a graph in coordinate form. Make the indexing 1-based */
void VSP_dump_graph
(
    VSP_state *X,
    VSP_state *Y,
    VSP_problem *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT col, row, i ;
    VSPLONG p ;

    /* Unpack graph */
    VSPINT *Ai = problem->Graph->Ai ;
    VSPLONG *Ap = problem->Graph->Ap ;
    VSPINT  *c = problem->Graph->c ;
    VSPINT  n  = problem->Graph->n  ;
    VSPLONG nz = problem->Graph->nz ;

    /* Open file */
    char path[1000];
    sprintf(path, "output/graph_%llu.txt", count);

    FILE *file = fopen(path, "w") ;

    /* Write n and nz */
    fprintf(file, "n: %lld\nnz: %lld\n", (long long) n, (long long) nz / 2) ;

    /* Write the cost values */
    fprintf(file, "\nCosts:") ;
    for(i = 0 ; i < n ; i++)
    {
        if(i + 1 < 10)
        {
            fprintf(file, "\nc(%lld):  %lld", (long long) i + 1, (long long) c[i]) ;
        }
        else
        {
            fprintf(file, "\nc(%lld): %lld", (long long) i + 1, (long long) c[i]) ;
        }
    }

    fprintf(file, "\n\nEdges:") ;
    for(col = 0 ; col < n ; col++)
    {
        for(p = Ap[col] ; p < Ap[col + 1] ; p++)
        {
            row = Ai[p] ;

            /* Print out only the bottom diagonal half of the matrix */
            if(row > col)
            {
                if(col + 1 < 10)
                {
                    fprintf(file, "\n%lld   %lld", (long long) col + 1, (long long) row + 1) ;
                }
                else
                {
                    fprintf(file, "\n%lld  %lld", (long long) col + 1, (long long) row + 1) ;
                }
            }

        }
    }

    if(X == NULL || Y == NULL)
    {
        fclose(file) ;
        return ;
    }

    /* Print x */
    VSPINT *x = X->point->int_part ;

    fprintf(file, "\n\nX:") ;
    for(i = 0 ; i < n ; i++)
    {
        if(x[i])
        {
            fprintf(file, "\n%lld", (long long) i + 1) ;
        }
    }

    /* Print y */
    VSPINT *y = Y->point->int_part ;
    fprintf(file, "\n\nY:") ;
    for(i = 0 ; i < n ; i++)
    {
        if(y[i])
        {
            fprintf(file, "\n%lld", (long long) i + 1) ;
        }
    }

    fclose(file) ;

    return ;
}


// Example usage:
// VSP_create_dot_ovi_file(problem, argv[file_num], "data/test_ovi/") ;
void VSP_create_dot_ovi_file
(
    VSP_problem *problem,
    char           *name,    /* From cmnd line. Ex: data / **** / bcspwr08.mtx */
    char         *outDir
)
{
    VSPLONG i ;
    VSPINT k ;
    k = 0 ; // To silence warning

    /* 1. Figure out file name */
    /* 1a. Set the prefix as data/ovi/*/
    char *prefix = outDir ;
    char file_name[500] ;
    file_name[0] = '\0' ;
    strcpy(file_name, prefix) ;

    /* 1b. Find last slash in name and copy after the slash, stopping 
           before the extension */
    VSPINT outDir_len = strlen(outDir);
    char *slash = strrchr(name, '/') ;
    for(i = 0 ; slash[i] != '.' ; i++ )
    {
        k = outDir_len + i ;
        file_name[k] = slash[i] ;
    }

    /* 1c. Add a .ovi extension */
    file_name[++k]   = '.' ;
    file_name[++k] = 'o' ;
    file_name[++k] = 'v' ;
    file_name[++k] = 'i' ;
    file_name[++k] = '\0' ;

    printf("\n%s", file_name) ;

    /* 2. Write the problem to a file */
    VSPINT *Ai = problem->Graph->Ai ;
    VSPLONG *Ap = problem->Graph->Ap ;
    VSPINT  n  = problem->Graph->n  ;
    VSPLONG nz = problem->Graph->nz ;

    /* 2a. Write n and nz */
    FILE *file = fopen(file_name, "w") ;
    fprintf(file, "n: %lld    nz: %lld", (long long) n, (long long) nz) ;

    /* 2b. Write Ap */
    fprintf(file, "\nAp:") ;
    for(i = 0 ; i <= n ; i++)
    {
        fprintf(file, "\n%lld", (long long) Ap[i]) ;
    }

    /* 2c. Write Ai */
    /* Leave one blank line between Ai and Ap */
    fprintf(file, "\n") ;
    fprintf(file, "\nAi:") ;
    for(i = 0 ; i < nz ; i++)
    {
        fprintf(file, "\n%lld", (long long) Ai[i]) ;
    }
    fprintf(file, "\n") ;

    fclose(file) ;

    return ;
}

// Write the partition to output/filename
// First, n is written as an int32
// Then, write a contiguous block of n int32s. The ith int32 represents the set
// in which node i is located.
// 0 means i is in X, i.e. x_i = 1  and  y_i = 0
// 1 means i is in Y, i.e. x_i = 0  and  y_i = 1
// 2 means i is in the separator, i.e. x_i = 1  and  y_i = 0
// 3 means i is in in the intersection of X and Y, i.e. x_i = 1  and  y_i = 1
void IO_Dump_Partition_int32
(
    VSPINT         *x,
    VSPINT         *y,
    VSPINT          n,
    char   *file_name  // Null terminated string. If NULL, then use default name
)
{
    // Open output/file_name in writing binary mode
    char path[1000] ;
    if(file_name == NULL) snprintf(path, 1000, "output/partition") ;
    else                  snprintf(path, 1000, "output/%s", file_name) ;
    FILE *fp = fopen(path, "wb") ;
    

    // Write n
    int32_t n_int32 = (int32_t) n ;
    fwrite(&n_int32, sizeof(int32_t), 1, fp) ;
    
    // Write the partition
    VSPINT i ;
    int32_t loc ;

    for(i = 0 ; i < n ; i++)
    {
        assert( IN_RANGE(x[i], 0, 2) && IN_RANGE(y[i], 0, 2) ) ;

        if     (x[i] == 1 && y[i] == 0) loc = 0 ;
        else if(x[i] == 0 && y[i] == 1) loc = 1 ;
        else if(x[i] == 0 && y[i] == 0) loc = 2 ;
        else if(x[i] == 1 && y[i] == 1) loc = 3 ;
        fwrite(&loc, sizeof(int32_t), 1, fp) ;
    }
    fclose(fp) ;
    return ;
}

void Print_Num_Flips(VSP_problem *P)
{
    static long long unsigned count = 0 ;
    count++ ;

    FILE *file ;
    if(count == 1)
    {
        file = fopen("output/Num_Flips.csv", "w") ;
        fprintf(file, "ProbName,n,nz,ProbLevel,NumFlips") ;
    }
    else file = fopen("output/Num_Flips.csv", "a") ;

    fprintf
    (
        file, "\n%s,%llu,%llu,%llu,%llu",
        P->name,
        CAST(P->Graph->n , long long unsigned),
        CAST(P->Graph->nz, long long unsigned),
        CAST(P->level, long long unsigned),
        CAST(P->stats->num_flips, long long unsigned)
    ) ;

    fclose(file) ;
    return ;
}


void cleanupGraphName(char *graph_name) ;

// Columns: ID  n  nz  sep  imb  time
void Write_Performance_Profile
(
    VSP_problem        *P,
    VSPINT       size_sep,
    VSPINT          wgt_x,
    VSPINT          wgt_y,
    VSPFLOAT         time
)
{
    static long long unsigned count = 0 ;
    count++ ;

    FILE *file ;
    char cpy_graph_name[1000] ;
    strcpy(cpy_graph_name, P->name) ;

    if(count == 1)
    {
        file = fopen("output/ppf.csv", "w") ;
        fprintf(file, "id,n,nz,sep,imb,time\n") ;
    }
    else file = fopen("output/ppf.csv", "a") ;

    /* Write data */
    fprintf(file, "%s,"  , cpy_graph_name) ;
    fprintf(file, "%lld,", (long long) P->Graph->n) ;
    fprintf(file, "%lld,", (long long) P->Graph->nz  ) ;
    fprintf(file, "%lld,", (long long) size_sep) ;
    fprintf(file, "%lld,", llabs(wgt_x - wgt_y)) ;
    fprintf(file, "%.8f" , time) ;
    fprintf(file, "\n") ;

    fclose(file) ;
    return ;
}
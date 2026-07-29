#include "core.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h>

/* Dump a graph as Gephi-ready CSVs: output/nodes_<count>.csv and
   output/edges_<count>.csv. Indexing is 1-based. If X and Y are
   non-NULL, node labels include the partition each node belongs to,
   based on the values of x[i] and y[i]. */
/* Dump a graph as Gephi-ready CSVs: output/nodes_<count>.csv and
   output/edges_<count>.csv. Indexing is 1-based. If X and Y are
   non-NULL, each node gets a low-cardinality Class column (for Gephi
   Partition coloring) plus a Label that carries the fractional value
   for the at-most-one fractional vertex per shore. The output/ dir
   must already exist. */
void dump_graph_gephi(Problem *problem, State *X, State *Y)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT col, row, i ;
    VSPLONG p ;

    /* Unpack graph */
    VSPINT *Ai = problem->graph->Ai ;
    VSPLONG *Ap = problem->graph->Ap ;
    VSPINT  *Ae = problem->graph->Ae ;
    VSPINT  *c  = problem->graph->c ;
    VSPINT  n   = problem->graph->n ;

    char path[1000] ;
    char label[100] ;
    char class[16] ;
    FILE *file ;

    /* ---------------- Nodes ---------------- */
    sprintf(path, "output/nodes_%llu.csv", count) ;

    file = fopen(path, "w") ;
    if(file == NULL) { perror(path) ; return ; }

    /* Hoist partition data: none of this depends on i */
    bool     labelled = (X != NULL && Y != NULL) ;
    VSPINT  *xint = NULL, *yint = NULL ;
    VSPINT   xind = -1, yind = -1 ;
    VSPFLOAT xfrac = 0, yfrac = 0 ;
    if(labelled)
    {
        xint  = X->point->int_part ;
        yint  = Y->point->int_part ;
        xind  = X->point->ind_frac ;
        yind  = Y->point->ind_frac ;
        xfrac = X->point->val_frac ;
        yfrac = Y->point->val_frac ;
    }

    fprintf(file, "Id,Label,Class,Cost\n") ;
    for(i = 0 ; i < n ; i++)
    {
        if(labelled)
        {
            bool xFrac = (i == xind) ;
            bool yFrac = (i == yind) ;

            if(xFrac && yFrac)
            {
                VSPFLOAT xval = xint[i] + xfrac ;
                VSPFLOAT yval = yint[i] + yfrac ;
                strcpy(class, "xy_frac") ;
                sprintf(label, "x_%.4f_cap_y_%.4f", xval, yval) ;
            }
            else if(xFrac)
            {
                VSPFLOAT xval = xint[i] + xfrac ;
                strcpy(class, "x_frac") ;
                sprintf(label, "x_%.4f", xval) ;
            }
            else if(yFrac)
            {
                VSPFLOAT yval = yint[i] + yfrac ;
                strcpy(class, "y_frac") ;
                sprintf(label, "y_%.4f", yval) ;
            }
            else if(xint[i] == 1 && yint[i] == 1)
            {
                strcpy(class, "xy") ;
                sprintf(label, "%lld", (long long) i + 1) ;
            }
            else if(xint[i] == 1)
            {
                strcpy(class, "x") ;
                sprintf(label, "%lld", (long long) i + 1) ;
            }
            else if(yint[i] == 1)
            {
                strcpy(class, "y") ;
                sprintf(label, "%lld", (long long) i + 1) ;
            }
            else
            {
                strcpy(class, "sep") ;
                sprintf(label, "%lld", (long long) i + 1) ;
            }
        }
        else
        {
            strcpy(class, "") ;
            sprintf(label, "%lld", (long long) i + 1) ;
        }

        fprintf(file, "%lld,%s,%s,%lld\n",
                (long long) i + 1,
                label,
                class,
                (long long) c[i]) ;
    }

    fclose(file) ;

    /* ---------------- Edges ---------------- */
    sprintf(path, "output/edges_%llu.csv", count) ;

    file = fopen(path, "w") ;
    if(file == NULL) { perror(path) ; return ; }

    fprintf(file, "Source,Target,Type,Weight\n") ;
    for(col = 0 ; col < n ; col++)
    {
        for(p = Ap[col] ; p < Ap[col + 1] ; p++)
        {
            row = Ai[p] ;

            /* Only the lower triangle, so each edge appears once */
            if(row > col)
            {
                fprintf(file, "%lld,%lld,Undirected,%lld\n",
                        (long long) col + 1,
                        (long long) row + 1,
                        (long long) Ae[p]) ;
            }
        }
    }

    fclose(file) ;

    return ;
}

/* Dump a graph in coordinate form. Make the indexing 1-based */
void dump_graph
(
    State *X,
    State *Y,
    Problem *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    VSPINT col, row, i ;
    VSPLONG p ;

    /* Unpack graph */
    VSPINT *Ai = problem->graph->Ai ;
    VSPLONG *Ap = problem->graph->Ap ;
    VSPINT  *c = problem->graph->c ;
    VSPINT  n  = problem->graph->n  ;
    VSPLONG nz = problem->graph->nz ;

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

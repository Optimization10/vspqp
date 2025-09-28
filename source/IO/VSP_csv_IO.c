/* This file is for writing numerical results to an existing .csv file        */
#include "shared_api.h"
#include <string.h>
#include <unistd.h>

#define BUFFER_SIZE 100000
#define MAX_LEN_GRAPH_NAME 1000

static void getColTitle
(
    char           *titleGiven,
    char         *titleToWrite
) ;

static void cleanupGraphName(char *graph_name) ;

/*----------------------------------------------------------------------------*/
/*   Write to Numerical_Results.csv                                           */
/*----------------------------------------------------------------------------*/

void VSP_write_to_csv
(
    char                *graph_name,
    char                *col_title_given, /* If empty string, a title is
                                             generated from the options       */
    VSP_solution_stats  *solStats,
    VSP_options         *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    char col_title_actual[500] = {'\0'} ;
    getColTitle(col_title_given, col_title_actual) ;

    char cpy_graph_name[MAX_LEN_GRAPH_NAME] ;
    strcpy(cpy_graph_name, graph_name) ;
    cleanupGraphName(cpy_graph_name) ;

    FILE *file ;
    if(count == 1)
    {
        file = fopen("output/Numerical_Results.csv", "w") ;

        /* Write column titles */
        fprintf(file, "Graph, " ) ;
        fprintf(file, "Vert, ") ;
        fprintf(file, "Edges, ") ;
        fprintf(file, "%s,", col_title_actual) ;
        fprintf(file, "Time\n") ;
    }
    else
    {
        file = fopen("output/Numerical_Results.csv", "a") ;
    }

    /* Write data */
    fprintf(file, "%s, ", cpy_graph_name) ;
    fprintf(file, "%i, ", solStats->num_nodes_total) ;
    fprintf(file, "%i, ", solStats->num_edges_total) ;
    fprintf(file, "%i, ", solStats->num_nodes_sep) ;
    fprintf(file, "%.5f\n", solStats->time) ;


    fclose(file) ;
    return ;
}




void VSP_write_to_csv_2
(
    char                *graph_name,
    char                *col_title_given, /* If empty string, a title is
                                             generated from the options       */
    VSP_solution_stats  *solStats,
    VSP_problem          *problem,
    VSP_options          *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    char col_title_actual[500] = {'\0'} ;
    getColTitle(col_title_given, col_title_actual) ;

    char cpy_graph_name[MAX_LEN_GRAPH_NAME] ;
    strcpy(cpy_graph_name, graph_name) ;
    cleanupGraphName(cpy_graph_name) ;

    FILE *file ;
    if(count == 1)
    {
        file = fopen("output/Numerical_Results.csv", "w") ;

        /* Write column titles */
        fprintf(file, "Graph," ) ;
        fprintf(file, "Vert,") ;
        fprintf(file, "Edges,") ;
        fprintf(file, "%s,", col_title_actual) ;
        fprintf(file, "Imbal,") ;
        fprintf(file, "Time\n") ;
    }
    else
    {
        file = fopen("output/Numerical_Results.csv", "a") ;
    }

    /* Write data */
    fprintf(file, "%s,", cpy_graph_name) ;
    fprintf(file, "%i,", solStats->num_nodes_total) ;
    fprintf(file, "%i,", solStats->num_edges_total) ;
    fprintf(file, "%i,", solStats->num_nodes_sep) ;
    fprintf(file, "%i,", abs(solStats->num_nodes_x - solStats->num_nodes_y)) ;
    fprintf(file, "%.6f\n", solStats->time) ;

    fclose(file) ;
    return ;
}














static void getColTitle
(
    char           *titleGiven,
    char         *titleToWrite
)
{
    if(*titleGiven == '\0')
    {
        strcpy(titleToWrite, "Ovi ") ;
    }
    else
    {
        if( strchr(titleGiven, ',') != NULL)
        {
            printf("\nError: No comma allowed in column title.\n") ;
            exit(1) ;
        }
        strcpy(titleToWrite, titleGiven) ;
    }
    return ;
}

/* The cleanup is done in-place */
static void cleanupGraphName(char *graph_name)
{
    /* Remove directories */
    char *last_slash = strrchr(graph_name, '/');
    if (last_slash != NULL)
    {
        memmove(graph_name, last_slash + 1, strlen(last_slash + 1) + 1);
    }

    /* Remove extension */
    char *last_dot = strrchr(graph_name, '.');
    if (last_dot != NULL)
    {
        *last_dot = '\0';
    }

    return ;
}
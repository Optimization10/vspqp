/* This file is for writing numerical results to an existing .csv file        */
#include "shared_api.h"
#include <string.h>
#include <unistd.h>

#define BUFFER_SIZE 100000
#define MAX_LEN_GRAPH_NAME 1000

void cleanupGraphName(char *graph_name) ;

/*----------------------------------------------------------------------------*/
/*   Write to Numerical_Results.csv                                           */
/*----------------------------------------------------------------------------*/









   





          



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
void cleanupGraphName(char *graph_name)
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
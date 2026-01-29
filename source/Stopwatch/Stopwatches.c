#include "Stopwatch.h"
#include "shared_api.h"
#include <stdlib.h>
#include "Stopwatches.h"
#include <string.h>

typedef struct Stopwatches
{
    size_t        numWatches ;
    size_t        maxNumWatches ;
    Stopwatch   **watches ;
    VSP_problem  *problem ; // Associated with this collection of stopwatches
} Stopwatches ;

Stopwatches *Stopwatches_Construct(size_t maxNumWatches, VSP_problem *problem)
{
    Stopwatches *SWS = calloc(1, sizeof(*SWS)) ;

    // SWS->numWatches    = 0 ;
    // SWS->maxNumWatches = maxNumWatches ;
    // SWS->watches       = calloc(maxNumWatches, sizeof(*SWS->watches)) ;
    // SWS->problem       = problem ;

    // // Make all watches NULL to be safe
    // for(size_t i = 0 ; i < maxNumWatches ; i++)
    // {
    //     SWS->watches[i] = NULL ;
    // }

    return SWS ;
}

void Stopwatches_Destroy(Stopwatches *SWS)
{
    for(size_t i = 0 ; i < SWS->numWatches ; i++)
    {
        Stopwatch_Destroy(SWS->watches[i]) ;
    }
    free(SWS->watches) ;
    free(SWS) ;

    // Do NOT destroy problem
}



void Stopwatches_Start(Stopwatches *SWS, char *name)
{
    // bool found = false ;
    // size_t whichWatch, i ;

    // // Search for a stopwatch with the same name
    // for(i = 0 ; i < SWS->numWatches ; i++)
    // {
    //     Stopwatch *SW = SWS->watches[i] ;
    //     if(strcmp( Stopwatch_GetName(SW), name) == 0)
    //     {
    //         found = true ;
    //         whichWatch = i ;
    //         break ;
    //     }
    // }

    // // Create a new stopwatch if it doesn't exist
    // if(!found)
    // {
    //     // Do not go over the maxNumWatches limit
    //     if(SWS->numWatches == SWS->maxNumWatches)
    //     {
    //         printf("\nERROR: Creating too many stopwatches") ;
    //         return ;
    //     }

    //     whichWatch = SWS->numWatches ;
    //     SWS->watches[SWS->numWatches++] = Stopwatch_Construct(name) ;
    // }

    // // Start the stopwatch
    // Stopwatch_Start(SWS->watches[whichWatch]) ;
    // return ;
}


void Stopwatches_Stop(Stopwatches *SWS, char *name)
{
    // bool found = false ;
    // size_t whichWatch, i ;

    // // Search for a stopwatch with the same name
    // for(i = 0 ; i < SWS->numWatches ; i++)
    // {
    //     Stopwatch *SW = SWS->watches[i] ;
    //     if(strcmp( Stopwatch_GetName(SW), name) == 0)
    //     {
    //         found = true ;
    //         whichWatch = i ;
    //         break ;
    //     }
    // }

    // if(!found)
    // {
    //     printf("\nERROR: Trying to end stopwatch %s which was not created",name);
    //     return ;
    // }

    // // Stop the stopwatch
    // Stopwatch_Stop(SWS->watches[whichWatch]) ;
}



void Stopwatches_Print(Stopwatches *SWS)
{
    // static bool firstRun = true;
    // FILE *file = fopen("output/stopwatches.csv", firstRun ? "w" : "a");
    // if (!file)
    // {
    //     printf("Error: Could not open stopwatches.csv\n");
    //     return;
    // }

    // if (firstRun)
    // {
    //     fprintf(file, "ProbName,n,nz,ProbLevel,SWname,SWtime\n");
    //     firstRun = false;
    // }

    // for (size_t i = 0; i < SWS->numWatches; i++)
    // {
    //     Stopwatch *sw = SWS->watches[i];
    //     fprintf
    //     (
    //         file, "%s,%lld,%lld,%lld,%s,%.6Lf\n",
    //         SWS->problem->name,
    //         CAST(SWS->problem->Graph->n, long long),
    //         CAST(SWS->problem->Graph->nz, long long),
    //         CAST(SWS->problem->level, long long),
    //         Stopwatch_GetName(sw),
    //         Stopwatch_GetAccumTime(sw)
    //     ) ;
    // }

    // fclose(file);
}
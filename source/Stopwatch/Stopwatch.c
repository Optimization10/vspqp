#include"shared_api.h"
#include<string.h>
#include"Stopwatch.h"


typedef struct Stopwatch
{
    char           *name ;
    struct timespec start ;
    struct timespec end ;
    bool            inProgress ;
    long double     accumTime ;
    long long int   counter ; // Number of times the sw has been start/stopped
} Stopwatch ;

Stopwatch *Stopwatch_Construct(char *SWname)
{
    // Create stopwatch with everything set to 0
    Stopwatch *SW = calloc(1, sizeof(*SW)) ;

    // Copy the name
    SW->name = malloc(strlen(SWname)+1) ;
    strcpy(SW->name, SWname) ;

    // Set everything to zero just to be safe
    SW->inProgress = false ;
    SW->accumTime  = 0 ;
    SW->counter    = 0 ;

    return SW ;
}

void Stopwatch_Destroy(Stopwatch *SW)
{
    free(SW->name) ;
    free(SW) ;
}


void Stopwatch_Start(Stopwatch *SW)
{
    // Input validation
    if(SW->inProgress)
    {
        printf("\nError: starting stopwatch when in progress") ;
        return ;
    }

    SW->inProgress = true ;
    SW->counter++ ;
    clock_gettime(CLOCK_MONOTONIC, &(SW->start)) ;
    return ;
}

void Stopwatch_Stop(Stopwatch *SW)
{
    clock_gettime(CLOCK_MONOTONIC, &(SW->end)) ;

    // Input validation
    if(!SW->inProgress)
    {
        printf("\nError: stopping an inactive stopwatch") ;
        return ;
    }

    SW->inProgress = false ;
    long double secDiff, nanoSecDiff ;

    secDiff        = CAST(SW->end.tv_sec - SW->start.tv_sec, long double) ;
    nanoSecDiff    = CAST(SW->end.tv_nsec - SW->start.tv_nsec, long double) ;
    SW->accumTime += secDiff + nanoSecDiff * CAST(1E-9, long double) ;
}

const char *Stopwatch_GetName(Stopwatch *SW)
{
    return SW->name ;
}

bool Stopwatch_GetInProgress(Stopwatch *SW)
{
    return SW->inProgress;
}

long double Stopwatch_GetAccumTime(Stopwatch *SW)
{
    return SW->accumTime;
}

long long int Stopwatch_GetCounter(Stopwatch *SW)
{
    return SW->counter;
}
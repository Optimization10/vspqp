#ifndef STOPWATCHES_H
#define STOPWATCHES_H

#include "Stopwatch.h"
#include "shared_api.h"

typedef struct Stopwatches Stopwatches ;

Stopwatches *Stopwatches_Construct(size_t maxNumWatches, VSP_problem *problem) ;
void Stopwatches_Destroy(Stopwatches *SWS) ;
void Stopwatches_Start(Stopwatches *SWS, char *name) ;
void Stopwatches_Stop(Stopwatches *SWS, char *name) ;
void Stopwatches_Print(Stopwatches *SWS) ;

#endif
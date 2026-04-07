#ifndef STOPWATCH_H
#define STOPWATCH_H
#include "core.h"
#include<stdbool.h>

typedef struct Stopwatch Stopwatch ;

Stopwatch *Stopwatch_Construct(char *SWname) ;
void Stopwatch_Destroy(Stopwatch *SW) ;
void Stopwatch_Start(Stopwatch *SW) ;
void Stopwatch_Stop(Stopwatch *SW) ;
const char *Stopwatch_GetName(Stopwatch *SW) ;
bool Stopwatch_GetInProgress(Stopwatch *SW) ;
long double Stopwatch_GetAccumTime(Stopwatch *SW) ;
long long int Stopwatch_GetCounter(Stopwatch *SW) ;
#endif 
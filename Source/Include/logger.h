#ifndef LOGGER_H
#define LOGGER_H

#include "dynamic_stack.h"

struct Logger
{
    bool enabled ;
    DynamicStack *initRecords ;
    FILE *log ;
} ;

typedef struct Logger_Record Logger_Record ;
typedef struct Logger Logger ;

extern Logger LOG_G ;


void LoggerInit(bool enable) ;
void LoggerDestroy() ;
void LoggerStart
(
    State *X, State *Y, Problem *problem,
    const char *format, ...
) ;
void LoggerStop(State *X, State *Y, Problem *problem) ;

#endif
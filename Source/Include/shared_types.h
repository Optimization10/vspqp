#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include<stdbool.h>
#include "core.h"
#include"VSPQP.h"


typedef enum
{
    UNCOMPRESSED,
    COMPRESSED
} Compression_Status ;


/* The algorithm always performs heavy edge matching and brotherly matching.
   This struct records the user's decision about when to enable community
   matching */



void LoggerInit(bool enable);
void LoggerShutdown(void);
void LoggerStart(VSP_state *Xinit, VSP_state *Yinit, VSP_problem *problem, const char *format, ...);
void LoggerEnd(VSP_state *Xfin, VSP_state *Yfin, VSP_problem *problem, const char *format, ...);

#endif
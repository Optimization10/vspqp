#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include<stdbool.h>
#include "core.h"
#include"VSPQP.h"


typedef enum
{
 VERSION_DEFAULT,
 VERSION_1,
 VERSION_2,
 VERSION_3,
 VERSION_4,
 VERSION_5,
 VERSION_6,
 VERSION_7,
 VERSION_8
} VersionNum ;

typedef enum
{
    STRICT_INCREASE,
    WEAK_INCREASE,
    NO_CHANGE,
    WEAK_DECREASE,
    STRICT_DECREASE,
    UNSPECIFIED_EFFECT
} Obj_Effect ;

typedef enum
{
    N1_OPTIMAL,            // No feasible, strictly beneficial flip
    N1_WEAKLY_OPTIMAL,     // No feasible, strictly beneficial flip AND there is at least one neutral downflip
    N1_STRICTLY_OPTIMAL,   // No feasible, strictly beneficial flip AND there are zero neutral downflips
    NOT_N1_OPTIMAL         // At least one feasible, strictly beneficial flip
} VSP_N1_Optimality ;






















typedef enum
{
    FLIPS_DEC_LEVELS,        /* Specified using a levels struct */
    FLIPS_DEC_NO_COMM,        /* Do flips if the current level was obtained w/o community*/
    FLIPS_DEC_BOUNDARY    /* Do flips if problem->ratioMaxCostNodes > decBoundary */
} VSP_flips_dec_type ;

typedef struct
{
    VSP_flips_dec_type decType ;

    VSP_levels *levels ;
    VSPFLOAT decBoundary ;

} VSP_flips_dec ;





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
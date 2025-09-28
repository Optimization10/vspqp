#ifndef FLIPS_API_H
#define FLIPS_API_H

/* This file contains functions involving flips which are to be exposed to the
   outside */

#include "shared_api.h"
#include "ds_api.h"

typedef struct
{
    VSPINT  ind ;
    char    var ;
    VSPINT gain ;
    int     dir ;
} Flip ;


void Flips
(
    VSP_state               *X,
    VSP_state               *Y,
    VSP_problem       *problem,
    VSPINT        maxNumPasses,
    VersionNum         version
) ;

#endif
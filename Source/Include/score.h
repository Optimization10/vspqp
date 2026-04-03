#ifndef SCORE_H
#define SCORE_H

#include "shared_api.h"


typedef struct
{
    VSPFLOAT obj ;
    VSPINT imbal ;
} PointScore ;

void PointScore_SetObj(PointScore *Score, VSPFLOAT obj) ;
void PointScore_SetImb(PointScore *Score, VSPINT imb) ;

VSPINT PointScore_Compare(PointScore *Score1, PointScore *Score2) ;

void PointScore_Update
(
    PointScore *Score, State *X, State *Y, Problem *problem
) ;

bool PointScore_IsStrictlyBetter(PointScore *New, PointScore *Base) ;

#endif
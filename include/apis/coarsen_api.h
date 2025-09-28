#ifndef COARSENING_API_H
#define COARSENING_API_H

#include "shared_api.h"


VSP_problem* Problem_Coarsen
(
    VSP_problem *P,
    VSP_options *options
) ;

VSP_problem* VSP_coarsen_all_the_way_down
(
    VSP_problem* problem,
    VSP_options* options
) ;

void VSP_refine
(
    VSP_problem* finegraph,
    VSP_state* fineX,
    VSP_state* fineY,
    VSP_problem* coarsegraph,
    VSP_state* coarseX,
    VSP_state* coarseY 
) ;

void Uncoarsen_Point
(
    VSP_problem    *fine_graph,
    VSP_state          *fine_X,
    VSP_state          *fine_Y,
    VSP_state        *coarse_X,
    VSP_state        *coarse_Y
) ;

void VSP_DecideComm
(
    VSP_problem       *problem,
    VSP_options       *options
) ;

#endif
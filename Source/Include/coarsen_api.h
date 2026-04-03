#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "shared_api.h"

void Graph_Simple_Match
(
    VSP_Graph              *Graph,
    VSPINT                    lvl,
    VSPINT              *matching,
    VSPINT              *matchmap,
    VSPINT           *invmatchmap,
    VSP_match_type     *matchtype,
    VSPINT                    *cn,
    VSP_options          *options
) ;

VSP_problem* Problem_Coarsen
(
    VSP_problem *P,
    VSP_options *options
) ;

void Uncoarsen_Point
(
    VSP_problem    *fine_graph,
    VSP_state          *fine_X,
    VSP_state          *fine_Y,
    VSP_state        *coarse_X,
    VSP_state        *coarse_Y
) ;

#endif
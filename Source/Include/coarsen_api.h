#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "core.h"
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

Problem* Problem_Coarsen
(
    Problem *P,
    VSP_options *options
) ;

void Uncoarsen_Point
(
    Problem    *fine_graph,
    State          *fine_X,
    State          *fine_Y,
    State        *coarse_X,
    State        *coarse_Y
) ;

Problem *Problem_Recursive_Coarsen
(
    Problem *probFinest,
    VSP_options    *options
) ;

#endif
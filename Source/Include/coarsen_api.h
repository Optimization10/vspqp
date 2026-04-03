#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "core.h"
#include "shared_api.h"

void Graph_Simple_Match
(
    Graph              *Graph,
    VSPINT                    lvl,
    VSPINT              *matching,
    VSPINT              *matchmap,
    VSPINT           *invmatchmap,
    match_type     *matchtype,
    VSPINT                    *cn,
    Options          *options
) ;

Problem* Problem_Coarsen
(
    Problem *P,
    Options *options
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
    Options    *options
) ;

#endif
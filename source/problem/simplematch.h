#ifndef SIMPLE_MATCH_H
#define SIMPLE_MATCH_H

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

#endif
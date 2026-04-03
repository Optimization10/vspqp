#ifndef DELETE_NODES_H
#define DELETE_NODES_H

#include "core.h"
#include "shared_api.h"


VSP_Graph *Compress_Graph
(
    VSP_Graph *GraphOrg,
    VSPINT        *keep,   // keep[i] = 1 if to keep i, 0 if to not keep i
    VSPINT    *newIndex    // Output
) ;

State *Compress_Binary_State
(
    VSP_Graph *compGraph,
    VSPINT     *newIndex,
    State      *Zorg
) ;

// Given a state for the compressed problem, return the corresponding state
// for the uncompressed problem where the vertices which were not kept are in
// the separator 
State *Uncompress_State
(
    VSP_Graph *compGraph,
    VSPINT     *newIndex,
    State      *Zcmp,
    VSPINT        *AiOrg,
    VSPLONG *ApOrg,
    VSPINT          nOrg
) ;


#endif
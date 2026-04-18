#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "core.h"


void Graph_Simple_Match
(
    Graph       *Graph,
    Matching        *M,
    VSPINT         lvl,
    Options   *options
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

int Check_Matrix
(
    VSPINT  *Ai,
    VSPLONG *Ap,
    VSPINT  *Ae,
    VSPINT   n,
    VSPLONG  nz,
    int      require_sorted,
    int      use_error_macro
) ;

#endif
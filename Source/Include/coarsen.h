#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "core.h"

#define LINEAR    0
#define QUADRATIC 1

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
    Problem   *probFinest,
    Options   *options,
    VSPINT   **sing_sol,
    VSPINT    *n_sings
) ;

void Set_Upper_Bound
(
    Problem *problem,
    Options *opts,
    int      interpolation
) ;

#endif
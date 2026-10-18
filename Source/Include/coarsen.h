#ifndef COARSENING_API_H
#define COARSENING_API_H
#include "core.h"

typedef enum
{
    LINEAR,
    QUADRATIC
} Interpolation ;

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

void Set_Upper_Bounds
(
    Problem      *problem,
    Options      *opts,
    VSPFLOAT      ratio_coarsest,
    Interpolation interpolation
) ;

#endif
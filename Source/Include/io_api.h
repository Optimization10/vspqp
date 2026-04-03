#ifndef IO_API_H
#define IO_API_H

#include "shared_api.h"


Problem* read_problem(char* file_name, Options* options) ;

/* Dump graph to a file in coordinate format */
void dump_graph
(
    State *X,
    State *Y,
    Problem *problem
) ;


void print_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    Problem *problem,
    char           *name,  /* name of the problem */
    VSPFLOAT        time
) ;


#endif
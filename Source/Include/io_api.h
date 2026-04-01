#ifndef IO_API_H
#define IO_API_H

#include "shared_api.h"


VSP_problem* VSP_read_problem(char* file_name, VSP_options* options) ;

/* Dump graph to a file in coordinate format */
void VSP_dump_graph
(
    VSP_state *X,
    VSP_state *Y,
    VSP_problem *problem
) ;


void VSP_print_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    VSP_problem *problem,
    char           *name,  /* name of the problem */
    VSPFLOAT        time
) ;


#endif
#ifndef IO_H
#define IO_H
#include "core.h"

void dump_graph_gephi
(
    Problem *problem
) ;

void dump_graph
(
    State *X,
    State *Y,
    Problem *problem
) ;

void IO_Dump_Partition_int32
(
    VSPINT         *x,
    VSPINT         *y,
    VSPINT          n,
    char   *file_name  // Null terminated string. If NULL, then use default name
) ;

#endif
#ifndef GR_API_H
#define GR_API_H

#include "shared_api.h"

typedef enum {
    div_by_two,
    use_minreduce
} VSP_reduction_type ;

VSPFLOAT VSP_minreduce(
    VSP_state         *X,
    VSP_state         *Y,
    VSP_options *options,
    VSP_problem *problem
) ;

#endif
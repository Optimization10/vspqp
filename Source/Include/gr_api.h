#ifndef GR_API_H
#define GR_API_H

#include "shared_api.h"

typedef enum {
    div_by_two,
    use_minreduce
} reduction_type ;

VSPFLOAT minreduce(
    State         *X,
    State         *Y,
    Options *options,
    Problem *problem
) ;

#endif
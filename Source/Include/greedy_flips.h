#ifndef GREEDY_FLIPS_H
#include "shared_api.h"
#include <stdbool.h>

void Greedy_Flips
(
    State     *X,
    State     *Y,
    Problem   *problem,
    VSPINT     maxNumPasses,
    bool       deep_search
) ;

#endif
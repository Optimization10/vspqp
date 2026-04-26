#ifndef GREEDY_FLIPS_H

#include "core.h"
#include <stdbool.h>

void Greedy_Flips
(
    State     *X,
    State     *Y,
    Problem   *problem,
    VSPINT     maxNumPasses,
    bool       deep_search,
    bool       set_frac
) ;

#endif
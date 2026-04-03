#ifndef FM_API_H
#define FM_API_H

#include "shared_api.h"
void FM_Metis
(
    State             *Z,
    State             *W,
    Problem     *problem,
    Options     *options,
    VSPINT        num_passes
) ;
#endif
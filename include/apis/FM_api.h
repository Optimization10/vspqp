#ifndef FM_API_H
#define FM_API_H

#include "shared_api.h"
void FM_Metis
(
    VSP_state             *Z,
    VSP_state             *W,
    VSP_problem     *problem,
    VSP_options     *options,
    VSPINT        num_passes
) ;
#endif
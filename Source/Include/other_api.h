#ifndef OTHER_API_H
#define OTHER_API_H

#include "shared_api.h"
void VSP_set_initial_guess(
    State *X,
    State *Y,
    Problem *problem,
    VSP_options *options
) ;

VSP_options* VSP_construct_options() ;


void VSP_preprocess_finest_problem
(
    Problem* problem,
    VSP_options* options
) ;


#endif
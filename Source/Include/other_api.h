#ifndef OTHER_API_H
#define OTHER_API_H

#include "shared_api.h"
void set_initial_guess(
    State *X,
    State *Y,
    Problem *problem,
    Options *options
) ;

Options* construct_options() ;


void preprocess_finest_problem
(
    Problem* problem,
    Options* options
) ;


#endif
#include "shared_api.h"
#include "other_api.h"
#include "coarsen_api.h"
#include "flips_api.h"

/*----------------------------------------------------------------------------*/
/*  Declarations                                                              */
/*----------------------------------------------------------------------------*/
void check_for_complete_graph(Problem *problem) ;

void ResolveLevelDecideAndDo
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options
) ;


void Optimize
(
    State                     *X,
    State                     *Y,
    Problem             *problem,
    Options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
) ;
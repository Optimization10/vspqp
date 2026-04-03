#include "shared_api.h"
#include "other_api.h"
#include "coarsen_api.h"
#include "flips_api.h"

/*----------------------------------------------------------------------------*/
/*  Declarations                                                              */
/*----------------------------------------------------------------------------*/
void VSP_check_for_complete_graph(Problem *problem) ;

void VSP_ResolveLevelDecideAndDo
(
    State             *X,
    State             *Y,
    Problem     *problem,
    VSP_options     *options
) ;


void Optimize
(
    State                     *X,
    State                     *Y,
    Problem             *problem,
    VSP_options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
) ;
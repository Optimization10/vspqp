#include "shared_api.h"
#include "other_api.h"
#include "coarsen_api.h"
#include "flips_api.h"

/*----------------------------------------------------------------------------*/
/*  Declarations                                                              */
/*----------------------------------------------------------------------------*/
void VSP_check_for_complete_graph(VSP_problem *problem) ;

void VSP_ResolveLevelDecideAndDo
(
    VSP_state             *X,
    VSP_state             *Y,
    VSP_problem     *problem,
    VSP_options     *options
) ;


void Optimize
(
    VSP_state                     *X,
    VSP_state                     *Y,
    VSP_problem             *problem,
    VSP_options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
) ;
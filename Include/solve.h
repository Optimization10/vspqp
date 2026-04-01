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

// whereToUseFMFlips:
// 0 - FMv2 disabled (or anything besides 1, 2, 3)
// 1 - FMv2 only during first optimization
// 2 - FMv2 during gamma reduction
// 3 - FMv2 in both locations
void Optimize
(
    VSP_state                     *X,
    VSP_state                     *Y,
    VSP_problem             *problem,
    VSP_options             *options,
    bool          useManyFlipsPoints,
    VSPINT         whereToUseFMFlips
) ;
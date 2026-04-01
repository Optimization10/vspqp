#include "core.h"
#include "shared_api.h"
#include "ds_api.h"
#include "coarsen_internal.h"



/* Modify options to reflect the decision */
void VSP_DecideComm
(
    VSP_problem       *problem,
    VSP_options       *options
)
{
    VSP_coarsen_dec *dec = options->coarsenDec ;

    if(dec == NULL)
    {
        options->do_community_match = (problem->level <= 1) ? false : true ;
        return ;
    }

    bool isMem, doComm ;
    switch(dec->decType)
    {
        case CRSN_DEC_useCommAtAllLevels:
            options->do_community_match = true ;
            break ;

        case CRSN_DEC_useCommAtNoLevels:
            options->do_community_match = false ;
            break ;

        case CRSN_DEC_useCommAtCertainLevels:
            isMem = VSP_levels_is_member(problem, dec->levels) ;
            options->do_community_match = isMem ;
            break ;

        case CRSN_DEC_useCommWithBoundary:
            doComm = ( problem->ratioMaxCostNodes < dec->decBounday ) ;
            options->do_community_match = doComm ;
            break ;

        default:
            printf("\n\nInvalid coarsening decision\n\n") ;
            exit(1) ;
    }
    return ;
}
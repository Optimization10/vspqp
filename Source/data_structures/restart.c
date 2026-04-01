#include "core.h"
#include "shared_api.h"
#include "ds_api.h" // IWYU pragma: keep





/* Just construct and set everything to empty */
VSP_restart_decision *VSP_construct_restart_decision()
{
    VSP_restart_decision *ret = malloc(sizeof(*ret)) ;

    /* Levels */
    ret->use_levels = false ;
    ret->levels     = NULL ;

    /* Rule */
    ret->use_rule             = false ;
    ret->restarts_done        = 0 ;
    ret->restarts_allowed     = 0 ;
    ret->rule_param           = 0 ;
    ret->finest_level_allowed = 0 ;

    return ret ;
}

void VSP_destruct_restart_decision(VSP_restart_decision *dec)
{
    if (dec->levels != NULL)
    {
        VSP_destruct_levels(dec->levels) ;
    }

    free(dec) ;

    return ;
}

/* If using the rule, and if it returns true, then it also increments
   dec->restarts_done. */
bool VSP_restart_decide
(
    VSP_restart_decision       *dec,
    VSPFLOAT                  objIn, /* If not using rule, this is ignored */
    VSPFLOAT                 objOut, /* If not using rule, this is ignored */
    VSP_problem            *problem  /* If using levels, this is ignored   */
)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        if(dec->use_levels + dec->use_rule != 1)
        {
            VSP_error("Error in restart decision", DBG_INFO) ;
        }

        if(dec->use_rule)
        {
            if(dec->restarts_done > dec->restarts_allowed)
            {
                VSP_error("Error in restart decision", DBG_INFO) ;
            }

            if(objOut < objIn)
            {
                VSP_error("Probable error. Final objective worse.", DBG_INFO) ;
            }
        }
        else
        {
            if(dec->levels == NULL)
            {
                VSP_error("Error in restart decision", DBG_INFO) ;
            }
        }
    #endif



    if (dec->use_levels)
    {
        return VSP_levels_is_member(problem, dec->levels) ;
    }
    else if (dec->use_rule)
    {
        if(dec->restarts_done >= dec->restarts_allowed
            || problem->level < dec->finest_level_allowed)
        {
            return false ;
        }

        /* Check if we should restart */
        VSPFLOAT rel_change = (objOut - objIn) / objIn ;
        if (rel_change <= dec->rule_param)
        {
            dec->restarts_done++ ;
            return true ;
        }
        else
        {
            return false ;
        }
    }

    return false ;
}

void VSP_restart_decide_reset(VSP_restart_decision *dec)
{
    dec->restarts_done = 0 ;
    return ;
}
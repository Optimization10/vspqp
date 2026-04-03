#include "core.h"
/*----------------------------------------------------------------------------*/
// states_manager module
// Purpose: In several parts of the program I want to start from some initial
// point (usually the best point so far), modify it, but have the option of
// reverting back. This module provides a container for the best and work
// points, as well as functions to manage these points
/*----------------------------------------------------------------------------*/
#include "shared_api.h"
#include "states_manager.h"

struct StatesMan
{
    State *best[2] ;
    State *work[2] ;
} ;




void StatesMan_Destroy(StatesMan *Man)
{
    Destroy_State(Man->work[0]) ;
    Destroy_State(Man->work[1]) ;
    free(Man) ;
    return ; 
}

// Deep copy best states into work states
void StatesMan_UpdateWorkStates(StatesMan *Man)
{
    copy_state(Man->best[0], Man->work[0]) ;
    copy_state(Man->best[1], Man->work[1]) ;
    return ;
}


// Deep copy work states into best states
void StatesMan_UpdateBestStates(StatesMan *Man)
{
    copy_state(Man->work[0], Man->best[0]) ;
    copy_state(Man->work[1], Man->best[1]) ;
    return ;
}

State *StatesMan_GetState
(
    StatesMan *SM, StateManStateType type, VSPINT which
)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(IN_RANGE(which, 0, 2)) ;

    switch (type)
    {
    case BEST_STATES :
        return SM->best[which] ;
    case WORK_STATES :
        return SM->work[which] ;
    default:
        return NULL ;
    }
}


StatesMan *StatesMan_Init(State *Xbest, State *Ybest)
{
    VSPINT n = Xbest->point->n ;
    StatesMan *ret = malloc(sizeof(StatesMan)) ;

    ret->best[0] = Xbest ;
    ret->best[1] = Ybest ;

    ret->work[0] = Construct_State(n) ;
    ret->work[1] = Construct_State(n) ;

    StatesMan_UpdateWorkStates(ret) ;
    
    return ret ;
}
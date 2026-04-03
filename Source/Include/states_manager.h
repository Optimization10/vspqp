#ifndef STATES_MANAGER_H
#define STATES_MANAGER_H
#include "shared_api.h"

typedef enum
{
    BEST_STATES,
    WORK_STATES
} StateManStateType ;

typedef struct StatesMan StatesMan ;
StatesMan *StatesMan_Init(State *Xbest, State *Ybest) ;
void StatesMan_Destroy(StatesMan *Man) ;
void StatesMan_UpdateWorkStates(StatesMan *Man) ;
void StatesMan_UpdateBestStates(StatesMan *Man) ;
State *StatesMan_GetState
(
    StatesMan *SM, StateManStateType type, VSPINT which
) ;
#endif
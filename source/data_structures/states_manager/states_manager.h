#ifndef STATES_MANAGER_H
#define STATES_MANAGER_H
#include "shared_api.h"

typedef enum
{
    BEST_STATES,
    WORK_STATES
} StateManStateType ;

typedef struct StatesMan StatesMan ;
StatesMan *StatesMan_Init(VSP_state *Xbest, VSP_state *Ybest) ;
void StatesMan_Destroy(StatesMan *Man) ;
void StatesMan_UpdateWorkStates(StatesMan *Man) ;
void StatesMan_UpdateBestStates(StatesMan *Man) ;
VSP_state *StatesMan_GetState
(
    StatesMan *SM, StateManStateType type, VSPINT which
) ;
#endif
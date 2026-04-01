#ifndef MEMORY_FLIPS_H
#define MEMORY_FLIPS_H

#include "shared_api.h"


typedef struct FlipsMemory FlipsMemory ;

FlipsMemory *FlipsMem_Construct(VSPINT capacity, FlipsCtx *FC) ;
void FlipsMem_Destroy(FlipsMemory *fm) ;
void FlipsMem_Update(FlipsMemory *FM, Flip *flip) ;
void FlipsMem_Rollback(FlipsMemory *FM) ;
void FlipsMem_Reset(FlipsMemory *FM) ;
VSPINT FlipsMem_GetBestAccumGain(FlipsMemory *FM) ;
#endif
#include "core.h"
// #include "source/data_structures/flips_memory/flips_mem.h"
#include "flips_context.h"
#include "flips_mem.h"

struct FlipsMemory
{
    // For i = 1,2,...numFlips, the ith flip is in memory[i-1] ; 
    Flip *memory ;
    VSPINT numFlips ;

    // Max number of flips we are capable of holding
    VSPINT capacity ;


    // Current accumulated gain and imbalance
    VSPINT currAccumGain ;
    VSPINT currImb ;

    // Point 0 is before any flip. Point 1 is after 1st flip. Point i is after
    // the ith flip for i = 1,2,3... Hence, when rolling back, we want to undo
    // flips nflip ,nflip-1 ... bestPoint+1, which are contained in
    // memory[nflip-1...bestPoint]
    VSPINT bestPoint ;
    VSPINT bestAccumGain ; // Starts off as 0
    VSPINT bestImb ;

    FlipsCtx *FC ;
} ;

FlipsMemory *FlipsMem_Construct(VSPINT capacity, FlipsCtx *FC)
{
    FlipsMemory *ret = calloc(1, sizeof(FlipsMemory)) ;
    ret->memory      = calloc(capacity, sizeof(Flip)) ;
    ret->capacity    = capacity ;
    ret->FC          = FC ;

    return ret ;
}

void FlipsMem_Destroy(FlipsMemory *fm)
{
    if (!fm) return;
    free(fm->memory);
    free(fm);
}

// Perform AFTER doing a flip
// ASSUMED: not going over capacity and not an empty flip
void FlipsMem_Update
(
    FlipsMemory    *FM,
    Flip         *flip
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef DEBUG
        if( (FM->numFlips == FM->capacity) || (flip->var == NO_VAR) )
        {   
            ERR("Problem in flips memory") ;
            return ;
        }
    #endif

    // Memory, current gain and imbalance
    FM->memory[FM->numFlips++] = *flip ;
    FM->currAccumGain += flip->gain ;
    FM->currImb = IMBALANCE( FlipsCtx_Get_Wgt(FM->FC, 0), FlipsCtx_Get_Wgt(FM->FC, 1) ) ;

    bool betterGain = FM->currAccumGain > FM->bestAccumGain ;
    bool sameGain   = FM->currAccumGain == FM->bestAccumGain ;

    if(betterGain)
    {
        // printf("\nBetterGain") ;
        FM->bestAccumGain = FM->currAccumGain ;
        FM->bestImb       = FM->currImb ;
        FM->bestPoint     = FM->numFlips ;
    }
    else if(sameGain)
    {
        if(FM->currImb <= FM->bestImb)
        {
            FM->bestAccumGain = FM->currAccumGain ;
            FM->bestImb       = FM->currImb ;
            FM->bestPoint     = FM->numFlips ;
        }
    }
}

// Roll back to the best point
void FlipsMem_Rollback(FlipsMemory *FM)
{
    // Undo flips in memory[nflip-1...bestPoint]
    VSPINT start = FM->numFlips - 1 ;
    VSPINT end   = FM->bestPoint ;
    for(VSPINT p = start ; p >= end ; p--)
    {
        FlipsCtx_Apply_Flip_Basic(FM->FC, &(FM->memory[p])) ;
    }
}


void FlipsMem_Reset(FlipsMemory *FM)
{
    FM->numFlips = 0;
    FM->currAccumGain = 0;
    FM->currImb = 0;
    FM->bestPoint = 0;
    FM->bestAccumGain = 0;
    FM->bestImb = 0;
}


VSPINT FlipsMem_GetBestAccumGain(FlipsMemory *FM)
{
    return FM->bestAccumGain ;
}
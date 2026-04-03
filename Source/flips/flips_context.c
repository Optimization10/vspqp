#include "core.h"
/*----------------------------------------------------------------------------*/
// flips_context module
// PURPOSE: Store and keep updated all the data pertaining to the point we
// are doing flips in. Also returns flips given various requirements
/*----------------------------------------------------------------------------*/


#include "shared_api.h"
#include "ds_api.h"
#include "flips_context.h"
#include "flips_common_api.h"

struct FlipsCtx
{
    VSPINT *vars [2] ;
    VSPINT *Avars[2] ;
    VSPINT  wgts [2] ;
    VSPINT *locked[2] ;
    Doubly_LLs_Sized *Lists[2][2] ;
    VSPINT Ts[2][2] ;
    VSPINT n ;
    Problem *problem ;

    bool randomize ; // randomize initialization of lists
} ;


/*----------------------------------------------------------------------------*/
// Definitions
/*----------------------------------------------------------------------------*/
// Ignore locked vertices.
// If randomize is true, then build the lists in a randomized manner by going
// through the vertices in a random order. Otherwise, go through the vertices
// in numerical order
//
// locked - If null, then ignored. One or more of its entries can also be NULL,
//          in which case that particular array is ignored




void FlipsCtx_Check
(
    FlipsCtx *FC
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT v, zi, T, Tdbg, k ;

    VSPINT *vars[2];
    VSPINT *Avars[2];
    VSPINT wgts[2];
    VSPINT *locked[2];
    VSPINT Ts[2][2];
    Doubly_LLs_Sized *Lists[2][2];

    for (v = 0 ; v < 2 ; v++)
    {
        vars[v]   = FC->vars[v];
        Avars[v]  = FC->Avars[v];
        wgts[v]   = FC->wgts[v];
        locked[v] = FC->locked[v];

        for (zi = 0; zi < 2; ++zi)
        {
            Ts[v][zi]    = FC->Ts[v][zi];
            Lists[v][zi] = FC->Lists[v][zi];
        }
    }

    Problem *problem          = FC->problem ;
    VSPINT n = problem->Graph->n ;

    /*------------------------------------------------------------------------*/
    // Check vars, Avars, weights, and feasibility by placing into states
    /*------------------------------------------------------------------------*/
    VSP_point
    ptX = {.int_part = vars[Z_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n},
    ptY = {.int_part = vars[W_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n} ;

    VSP_Az
    structAx = {.int_part   = Avars[Z_VAR] , .ind_frac_z = EMPTY , .n = n   ,
                .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   },
    structAy = {.int_part   = Avars[W_VAR] , .ind_frac_z = EMPTY , .n = n   ,
                .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   };

    VSP_state
    X = {.point = &ptX, .Az = &structAx, .weight = wgts[Z_VAR]},
    Y = {.point = &ptY, .Az = &structAy, .weight = wgts[W_VAR]} ;
    VSP_check_states(&X, &Y, problem) ;

    /*------------------------------------------------------------------------*/
    // Check linked lists
    /*------------------------------------------------------------------------*/
    for(v = 0 ; v < 2 ; v++)
    {
        Doubly_LLs_Sized_Check
        (
            Lists[v][0], Lists[v][1] , vars[v], Avars[OTHER_VAR(v)], locked[v]
        ) ;
    }

    /*------------------------------------------------------------------------*/
    // Check Ts
    /*------------------------------------------------------------------------*/
    VSPINT ps[4], cs[4], wit[4] ;
    Compute_ps_And_cs_Dbg
    (
        vars[Z_VAR], Avars[Z_VAR], wgts[Z_VAR], locked[Z_VAR],
        vars[W_VAR], Avars[W_VAR], wgts[W_VAR], locked[W_VAR],
        problem, ps, cs, wit
    ) ;

    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            T = Ts[v][zi] ;
            k = GET_ENCODING(v, zi) ;
            Tdbg = ps[k] ;
            if(T != Tdbg)
            {
                ERR("T incorrect") ;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Check locked arrays
    /*------------------------------------------------------------------------*/
    // Check that they are 0, 1
    for(v = 0 ; v < 2 ; v++) Check_Array_Is_Binary(locked[v], n) ;

    return ;
}

void FlipsCtx_Initialize_Lists_And_Ts(FlipsCtx *FC)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i, v, j,*perm, zi, Awi, encoding, dir,
    n = FC->problem->Graph->n ;

    bool randomize = FC->randomize ;

    VSPINT *vars[2];
    VSPINT *Avars[2];
    VSPINT wgts[2];
    VSPINT *locked[2] ;

    for (v = 0 ; v < 2 ; v++)
    {
        vars[v]   = FC->vars[v];
        Avars[v]  = FC->Avars[v];
        wgts[v]   = FC->wgts[v];
        locked[v] = FC->locked[v];
    }



    /*------------------------------------------------------------------------*/
    // Reset lists
    /*------------------------------------------------------------------------*/
    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            Doubly_LLs_Sized_Reset(FC->Lists[v][zi]) ;
        }
    }

    if(randomize)
    {
        perm = malloc(n * sizeof(*perm)) ;
        Get_Permutation(perm, n) ;
    }

    /*------------------------------------------------------------------------*/
    // Initialize
    /*------------------------------------------------------------------------*/
    // TODO: If no heap, use D+2 where D is max degree
    for(v = 0 ; v < 2 ; v++)
    {
        FC->Ts[v][0] = n + 1 ;
        FC->Ts[v][1] = -1 ;
    }

    for(j = 0 ; j < n ; j++)
    {
        i = randomize ? perm[j] : j ;
        for(v = 0 ; v < 2 ; v++)
        {
            // Skip locked vertices
            if( !UNLOCKED(locked[v], i) ) continue ;

            zi  = vars[v][i] ;
            Awi = Avars[GET_OTHER_VAR(v)][i] ;
            dir = DIR_FLIP(zi) ;

            // Insert into ll and update T
            Doubly_LLs_Sized_Insert(FC->Lists[v][zi], Awi, i);
            if(dir * Awi < dir * FC->Ts[v][zi] ) FC->Ts[v][zi] = Awi ;
        }
    }

    /*------------------------------------------------------------------------*/
    // Cleanup
    /*------------------------------------------------------------------------*/
    if(randomize)
    {
        free(perm) ;
    }
    return ;
}


FlipsCtx *FlipsCtx_Init
(
    VSP_state   *states[2],
    bool         randomize,
    Problem   *problem
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT v, zi,
    n = problem->Graph->n ;

    FlipsCtx *FC = malloc( sizeof(FlipsCtx) ) ;
    if(!FC) return NULL ;

    VSP_state *Z = states[0], *W = states[1] ;

    // Vars, Avars, wgts
    FC->vars [0] = Z->point->int_part ;
    FC->Avars[0] = Z->Az->int_part ;
    FC->wgts [0] = Z->weight ;

    FC->vars [1] = W->point->int_part ;
    FC->Avars[1] = W->Az->int_part ;
    FC->wgts [1] = W->weight ;

    FC->randomize = randomize ;

    for(v = 0 ; v < 2 ; v++) FC->locked[v] = calloc(n, sizeof(VSPINT)) ;

    // Allocate space for lists
    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            FC->Lists[v][zi] = Doubly_LLs_Sized_Init(n) ;
        }
    }

    // Problem and n
    FC->problem = problem ;
    FC->n = n ;

    // Initialize lists and Ts
    FlipsCtx_Initialize_Lists_And_Ts(FC) ;

    IFDBG(FlipsCtx_Check(FC)) ;
    return FC ;
}

void FlipsCtx_Destroy(FlipsCtx *FC)
{
    if(!FC) return ;

    VSPINT v, zi ;
    for(v = 0 ; v < 2 ; v++) free(FC->locked[v]) ;
    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            Doubly_LLs_Sized_Destroy(FC->Lists[v][zi]) ;
        }
    }
    free(FC) ;
    return ;
}

// Re-compute the lists and Ts. Reset locked array
void FlipsCtx_Reset
(
    FlipsCtx *FC
)
{
    for(VSPINT i = 0 ; i < FC->n ; i++)
    {
        FC->locked[0][i] = 0 ;
        FC->locked[1][i] = 0 ;
    }

    FlipsCtx_Initialize_Lists_And_Ts(FC) ;

    IFDBG(FlipsCtx_Check(FC)) ;
    return ;
}

static bool FlipsCtx_TExists(FlipsCtx *FC, VSPINT T, VSPINT zi)
{
    IFDBG(FlipsCtx_Check(FC)) ;
    if(zi) return (T != -1) ;
    else   return (T != FC->n + 1) ;
}

// Return feasibility status. Also returns false if this flip DNE
// Get the gain of the flip in the entry of the optimal linked list for the pair
// v, zi. If the flip is infeasible or does not exist for some reason, then
// feasible is set to false and gain may be garbage
static bool FlipsCtx_GetFeasStatusAndGain
(
    FlipsCtx       *FC,
    VSPINT           v,
    VSPINT          zi,
    VSPFLOAT     gamma,
    VSPFLOAT     *gain  // if feasible is false, this is garbage
)
{
    IFDBG(FlipsCtx_Check(FC)) ;
    VSPINT i, ci,
    T = FC->Ts[v][zi] ;

    // Determine if there is any flip corresponding to v, zi whatsoever
    if( !FlipsCtx_TExists(FC, T, zi) )
    {
        return false ;
    }

    // Determine feasibility
    Doubly_LLs_Sized_See_Entry(FC->Lists[v][zi], T, &i) ;
    ci = FC->problem->Graph->c[i] ;
    if(zi) // downflip
    {
        if(FC->wgts[v] == ci) return false ;
    }
    else // upflip
    {
        if(FC->wgts[v] + ci > FC->problem->u) return false ;
    }
    
    *gain = GAIN_FLIP(1 - 2 * zi, ci, gamma, T) ;
    return true ;
}


static void FlipsCtx_InitializeRuledIn(bool ruledIn[2][2])
{
    for (VSPINT v = 0 ; v < 2 ; v++)
    {
        for (VSPINT zi = 0 ; zi < 2 ; zi++)
        {
            ruledIn[v][zi] = false;
        }
    }
}


static void FlipsCtx_RuleInAFlip
(
    Flip_Selector_Mode           mode,
    VSPFLOAT                     gain,
    VSPINT                         zi,
    bool                ruledIn[2][2],
    VSPFLOAT              gains[2][2],
    VSPINT                          v
)
{
    if (mode == FEAS_FLIP)
    {
        ruledIn[v][zi] = true ;
        gains[v][zi] = gain ;
    }
    else if (mode == FEAS_BEN_FLIP)
    {
        if (gain > 1E-8 || (VSP_equal(gain, 0) && zi == 1))
        {
            ruledIn[v][zi] = true;
            gains[v][zi] = gain;
        }
    }
}

// Return true if we ruled in ANY flip whatsoever
static bool FlipsCtx_RuleInFlips
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars,
    Flip_Selector_Mode           mode,
    VSPFLOAT                    gamma,
    bool                ruledIn[2][2],
    VSPFLOAT              gains[2][2],
    VSPFLOAT                 *maxGain
)
{
    VSPINT v, zi ;
    VSPFLOAT gain = 0 ;
    bool feasible, correctVar, ret = false ;

    *maxGain = -(FC->problem->Graph->n * FC->problem->Graph->max_vtx_cost + FC->problem->Graph->max_vtx_cost) ;

    // First, rule out all flips
    FlipsCtx_InitializeRuledIn(ruledIn) ;

    // Rule in v, zi
    for(v = 0 ; v < 2 ; v++)
    {
        // Do not rule in any flips of the wrong variable
        correctVar = (whichVars == BOTH_VARS || whichVars == v) ;
        if(!correctVar) continue ;

        for(zi = 0 ; zi < 2 ; zi++)
        {
            feasible = FlipsCtx_GetFeasStatusAndGain(FC, v, zi, gamma, &gain);
            if(!feasible) continue ;

            // Depending on the mode, rule in a flip
            FlipsCtx_RuleInAFlip(mode, gain, zi, ruledIn, gains, v) ;
            if(!ret && ruledIn[v][zi]) ret = true ;
            
            // Update max feasible gain
            if(ruledIn[v][zi] && gain > *maxGain) *maxGain = gain ;
        }
    }
    return ret ;
}

// Return sum of weights
static void FlipsCtx_GetFlipWgts
(
    FlipsCtx             *FC,
    bool       ruledIn[2][2],
    VSPINT    flipWgts[2][2],
    VSPFLOAT     gains[2][2],
    VSPFLOAT         maxGain
)
{
    VSPINT v, zi ;
    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            if(ruledIn[v][zi] && gains[v][zi] == maxGain)
            {
                // Doubly_LLs_Sized *list = FC->Lists[v][zi] ;
                // flipWgts[v][zi] = Doubly_LLs_Sized_Get_Size(list, FC->Ts[v][zi]) ;
                flipWgts[v][zi] = 1 ;
            }
        }
    }
}


// Function to select a random entry based on weights
// Returns 1 if successful, 0 if all entries are 0
// Stores selected indices in v and col
static VSPINT FlipsCtx_Select_Weighted_Random(VSPINT wgts[2][2], VSPINT *v, VSPINT *zi)
{
    // Calculate sum of all entries
    VSPINT sum = 0 ;
    for (VSPINT i = 0 ; i < 2 ; i++)
    {
        for (VSPINT j = 0 ; j < 2 ; j++)
        {
            sum += wgts[i][j] ;
        }
    }
    
    // If sum is 0, no valid selection possible
    if (sum == 0)
    {
        return 0 ;
    }
    
    // Generate random number between 1 and sum
    VSPINT r = rand() % sum + 1 ;
    
    // Find the entry corresponding to the random number
    VSPINT currentSum = 0 ;
    for (VSPINT i = 0 ; i < 2 ; i++)
    {
        for (VSPINT j = 0 ; j < 2 ; j++)
        {
            currentSum += wgts[i][j] ;
            if (r <= currentSum)
            {
                *v  = i ;
                *zi = j ;
                return 1 ;
            }
        }
    }
    
    return 0 ; // Should never reach here
}

// whichVars = 0 --- var 0
// whichVars = 1 --- var 1
// whichVars = 2 --- both vars

// mode = 0 --- select the best feasible, even with negative gain 
// mode = 1 --- select the best beneficial flip. Return EMPTY if no
//              feas ben flip
Flip FlipsCtx_Select_Flip
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars,
    Flip_Selector_Mode           mode,
    VSPFLOAT                    gamma
)
{
    static long long unsigned count = 0 ; count++ ;
    TMP_ASSERT( IN_RANGE(whichVars, 0, 3) ) ;
    IFDBG(FlipsCtx_Check(FC)) ;

    VSPINT v, zi ;

    // Construct ruledIn
    bool ruledIn[2][2], any ;
    VSPFLOAT gains[2][2], maxGain ;

    any = FlipsCtx_RuleInFlips(FC, whichVars, mode, gamma, ruledIn, gains, &maxGain) ;
    if( !any )
    {
        // Return an empty flip
        return (Flip) {.ind = EMPTY, .var = NO_VAR} ;
    }

    // Construct weights
    VSPINT flipWgts[2][2] = {0} ;
    FlipsCtx_GetFlipWgts(FC, ruledIn, flipWgts, gains, maxGain) ;

    // Calculate probability
    VSPINT q = FlipsCtx_Select_Weighted_Random(flipWgts, &v, &zi) ;
    #ifdef VSP_DEBUG
    if(!q)
    {
        ERR("None selected") ;
    }

    if(!IN_RANGE(v, 0, 2) || !IN_RANGE(zi, 0, 2))
    {
        ERR("Out of bounds") ;
    }
    #endif

    VSPINT i,
    T = FC->Ts[v][zi] ;
    Doubly_LLs_Sized_See_Entry(FC->Lists[v][zi], T, &i) ;

    IFDBG(FlipsCtx_Check(FC)) ;
    return (Flip) {.ind = i, .var = v, .gain = maxGain} ;
}

/*----------------------------------------------------------------------------*/
// FM Flips v9 flip selector
/*----------------------------------------------------------------------------*/

static void FlipsCtx_GetFlipWgts_FMv9
(
    FlipsCtx             *FC,
    bool       ruledIn[2][2],
    VSPINT    flipWgts[2][2],
    VSPINT       gains[2][2],
    VSPINT           maxGain
)
{
    VSPINT v, zi ;
    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            if(ruledIn[v][zi] && gains[v][zi] == maxGain)
            {
                Doubly_LLs_Sized *list = FC->Lists[v][zi] ;
                flipWgts[v][zi] = Doubly_LLs_Sized_Get_Size(list, FC->Ts[v][zi]) ;
                // flipWgts[v][zi] = 1 ;
            }
        }
    }
}


// Return feasibility status. Also returns false if this flip DNE
// Get the gain of the flip in the entry of the optimal linked list for the pair
// v, zi. If the flip is infeasible or does not exist for some reason, then
// feasible is set to false and gain may be garbage
static bool FlipsCtx_GetFeasStatusAndGain_FMv9
(
    FlipsCtx       *FC,
    VSPINT           v,
    VSPINT          zi,
    VSPINT       *gain  // if feasible is false, this is garbage
)
{
    IFDBG(FlipsCtx_Check(FC)) ;
    VSPINT i, ci,
    T = FC->Ts[v][zi],
    gamma = FC->problem->Graph->max_vtx_cost ;

    // Determine if there is any flip corresponding to v, zi whatsoever
    if( !FlipsCtx_TExists(FC, T, zi) )
    {
        return false ;
    }

    // Determine feasibility
    Doubly_LLs_Sized_See_Entry(FC->Lists[v][zi], T, &i) ;
    ci = FC->problem->Graph->c[i] ;
    if(zi) // downflip
    {
        if(FC->wgts[v] == ci) return false ;
    }
    else // upflip
    {
        if(FC->wgts[v] + ci > FC->problem->u) return false ;
    }
    
    *gain = GAIN_FLIP(1 - 2 * zi, ci, gamma, T) ;
    return true ;
}

static void FlipsCtx_RuleInAFlip_FMv9
(
    VSPINT                       gain,
    VSPINT                         zi,
    bool                ruledIn[2][2],
    VSPINT                gains[2][2],
    VSPINT                          v
)
{
    if(zi == 1)
    {
        if(gain >= 0)
        {
            ruledIn[v][zi] = true ;
            gains[v][zi] = gain ;
        }
    }
    else
    {
        ruledIn[v][zi] = true ;
        gains[v][zi] = gain ;
    }
}

// Return true if we ruled in ANY flip whatsoever
static bool FlipsCtx_RuleInFlips_FMv9
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars,
    bool                ruledIn[2][2],
    VSPINT                gains[2][2],
    VSPINT                   *maxGain
)
{
    VSPINT v, zi ;
    VSPINT gain = 0, gamma = FC->problem->Graph->max_vtx_cost ;
    bool feasible, correctVar, ret = false ;

    *maxGain = -(FC->problem->Graph->n * FC->problem->Graph->max_vtx_cost + FC->problem->Graph->max_vtx_cost) ;

    // First, rule out all flips
    FlipsCtx_InitializeRuledIn(ruledIn) ;

    // Rule in v, zi
    for(v = 0 ; v < 2 ; v++)
    {
        // Do not rule in any flips of the wrong variable
        correctVar = (whichVars == BOTH_VARS || whichVars == v) ;
        if(!correctVar) continue ;

        for(zi = 0 ; zi < 2 ; zi++)
        {
            feasible = FlipsCtx_GetFeasStatusAndGain_FMv9(FC, v, zi, &gain);
            if(!feasible) continue ;

            // Depending on the mode, rule in a flip
            FlipsCtx_RuleInAFlip_FMv9(gain, zi, ruledIn, gains, v) ;
            if(!ret && ruledIn[v][zi]) ret = true ;
            
            // Update max feasible gain
            if(ruledIn[v][zi] && gain > *maxGain) *maxGain = gain ;
        }
    }
    return ret ;
}

// Selects the flip with best gain under the following constraints:
// feasible, unlocked downflip with gain >= 0
// feasible, unlocked upflip with any gain
Flip FlipsCtx_Select_Flip_FMv9
(
    FlipsCtx                      *FC,
    Flip_Selector_Vars      whichVars
)
{
    static long long unsigned count = 0 ; count++ ;
    TMP_ASSERT( IN_RANGE(whichVars, 0, 3) ) ;
    IFDBG(FlipsCtx_Check(FC)) ;

    VSPINT v, zi,
    gamma = FC->problem->Graph->max_vtx_cost ;

    // Construct ruledIn
    bool ruledIn[2][2], any ;
    VSPINT gains[2][2], maxGain ;

    any = FlipsCtx_RuleInFlips_FMv9(FC, whichVars, ruledIn, gains, &maxGain) ;
    if( !any )
    {
        // Return an empty flip
        return (Flip) {.ind = EMPTY, .var = NO_VAR} ;
    }

    // Construct weights
    VSPINT flipWgts[2][2] = {0} ;
    FlipsCtx_GetFlipWgts_FMv9(FC, ruledIn, flipWgts, gains, maxGain) ;

    // Calculate probability
    VSPINT q = FlipsCtx_Select_Weighted_Random(flipWgts, &v, &zi) ;
    #ifdef VSP_DEBUG
    if(!q)
    {
        ERR("None selected") ;
    }

    if(!IN_RANGE(v, 0, 2) || !IN_RANGE(zi, 0, 2))
    {
        ERR("Out of bounds") ;
    }
    #endif

    VSPINT i,
    T = FC->Ts[v][zi] ;
    Doubly_LLs_Sized_See_Entry(FC->Lists[v][zi], T, &i) ;

    IFDBG(FlipsCtx_Check(FC)) ;
    return (Flip) {.ind = i, .var = v, .gain = maxGain} ;
}


void FlipsCtx_Apply_Flip
(
    FlipsCtx      *FC,
    Flip        *flip,
    bool       lockIt
)
{
    SET_FUNCTION_COUNTER ;
    #ifdef VSP_DEBUG
        if (!FC || !FC->problem || !flip) {
            ERR("Invalid FlipsCtx, problem, or flip pointer");
            return;
        }
    #endif

    IFDBG(FlipsCtx_Check(FC));

    VSPLONG  p ;
    VSPINT j, T, wj, Azj, ziNew, zi, v;
    VSPINT i = flip->ind;
    VSPINT var = flip->var;
    VSPINT otherVar = OTHER_VAR(var);

    Problem *problem = FC->problem;
    VSPINT n = problem->Graph->n;
    VSPINT *z = FC->vars[var];
    VSPINT *Az = FC->Avars[var];
    VSPINT *w = FC->vars[otherVar];
    VSPINT *Aw = FC->Avars[otherVar];
    VSPINT *c = problem->Graph->c;
    VSPINT *Ai = problem->Graph->Ai;
    VSPLONG *Ap = problem->Graph->Ap;

    zi = z[i];
    VSPINT dir = DIR_FLIP(zi);
    VSPINT Awi = FC->Ts[var][zi];

    /* Remove i from its list and update T */
    Doubly_LLs_Sized_Delete(FC->Lists[var][zi], Awi, i);
    T = Awi;
    if (zi) { // downflip
        while (T > -1 && FC->Lists[var][zi]->sizes[T] == 0) T--;
    } else { // upflip
        while (T < n + 1 && FC->Lists[var][zi]->sizes[T] == 0) T++;
    }
    FC->Ts[var][zi] = T;

    /* Update zi, weight, and locked status */
    IFDBG({
        VSPINT ci = c[i];
        if (zi)
        {
            TMP_ASSERT(FC->wgts[var] > ci);
        } else
        {
            TMP_ASSERT(FC->wgts[var] + ci <= problem->u);
        }
    });
    z[i] += dir;
    FC->wgts[var] += dir * c[i];
    if (lockIt) FC->locked[var][i]++;

    /* Add i to new list if not locked */
    ziNew = zi + dir;
    if (!lockIt) {
        Doubly_LLs_Sized_Insert(FC->Lists[var][ziNew], Awi, i);
        T = FC->Ts[var][ziNew];
        if (dir * Awi > dir * T) FC->Ts[var][ziNew] = Awi;
    }

    /* Update neighbors */
    VSPLONG start = Ap[i];
    VSPLONG end = Ap[i + 1];
    p = start;
    j = i;
    bool iItselfFinished = false;
    do {
        if (iItselfFinished) {
            j = Ai[p];
            p++;
        }
        iItselfFinished = true;

        if (!FC->locked[otherVar][j]) {
            wj = w[j];
            Azj = Az[j];
            Doubly_LLs_Sized *lists = FC->Lists[otherVar][wj];
            Doubly_LLs_Sized_Delete(lists, Azj, j);
            Doubly_LLs_Sized_Insert(lists, Azj + dir, j);

            if (dir + 2 * wj == 1) {
                T = FC->Ts[otherVar][wj];
                if (lists->sizes[T] == 0) {
                    FC->Ts[otherVar][wj] += dir;
                }
            } else {
                if (Azj == FC->Ts[otherVar][wj]) {
                    FC->Ts[otherVar][wj] += dir;
                }
            }
        }
        Az[j] += dir;
    } while (p != end);

    IFDBG(FlipsCtx_Check(FC));
}

// Completely ignore everything except vars, Avars, and wgts
void FlipsCtx_Apply_Flip_Basic
(
    FlipsCtx      *FC,
    Flip        *flip
)
{
    #ifdef VSP_DEBUG
        Check_Weight(FC->vars[0], EMPTY, 0, FC->problem->Graph->c, FC->wgts[0], FC->problem->Graph->n) ;
        Check_Weight(FC->vars[1], EMPTY, 0, FC->problem->Graph->c, FC->wgts[1], FC->problem->Graph->n) ;
    #endif

    Do_Flip_Basic
    (
        FC->vars[flip->var],
        FC->Avars[flip->var],
        flip->ind,
        &(FC->wgts[flip->var]),
        FC->problem->Graph->c,
        FC->problem->Graph->Ai,
        FC->problem->Graph->Ap
    ) ;

    #ifdef VSP_DEBUG
        Check_Weight(FC->vars[0], EMPTY, 0, FC->problem->Graph->c, FC->wgts[0], FC->problem->Graph->n) ;
        Check_Weight(FC->vars[1], EMPTY, 0, FC->problem->Graph->c, FC->wgts[1], FC->problem->Graph->n) ;
    #endif
}


VSPINT **FlipsCtx_Get_Vars(FlipsCtx *FC)
{
    return FC->vars ;
}

VSPINT **FlipsCtx_Get_Avars(FlipsCtx *FC)
{
    return FC->Avars ;
}

VSPINT *FlipsCtx_Get_Wgts(FlipsCtx *FC)
{
    return FC->wgts ;
}

VSPINT **FlipsCtx_Get_Locked(FlipsCtx *FC)
{
    return FC->locked ;
}

Doubly_LLs_Sized *(*FlipsCtx_Get_Lists(FlipsCtx* FC))[2][2]
{
    return &FC->Lists;
}


VSPINT (*FlipsCtx_Get_Ts(FlipsCtx* FC))[2]
{
    return FC->Ts;
}

VSPINT FlipsCtx_Get_n(FlipsCtx* FC)
{
    return FC->n ;
}

Problem *FlipsCtx_Get_Problem(FlipsCtx* FC)
{
    return FC->problem ;
}

bool FlipsCtx_Is_Nonempty_Flip(Flip *flip)
{
    return (flip->ind != EMPTY) ;
}

// If this is an empty flip, return infeasible
// If u is empty, then problem->u is used
bool FlipsCtx_Is_Feasible_Flip(FlipsCtx *FC, Flip *flip, VSPINT u)
{
    SET_FUNCTION_COUNTER ;
    if(u == EMPTY) u = FC->problem->u ;

    if ( !FlipsCtx_Is_Nonempty_Flip(flip) ) return false ;

    VSPINT *z = FC->vars[flip->var] ;
    VSPINT *c = FC->problem->Graph->c ;
    VSPINT wgtVar = FC->wgts[flip->var] ;

    if(z[ flip->ind ])
    {
        return (wgtVar + c[flip->ind] <= u) ;
    }
    else
    {
        return (wgtVar > c[flip->ind]) ;
    }
}


VSPINT FlipsCtx_Get_Wgt
(
    FlipsCtx  *FC,
    VSPINT    var
)
{
    SET_FUNCTION_COUNTER ;
    TMP_ASSERT(var == 0 || var == 1) ;
    return FC->wgts[var] ;
}

void FlipsCtx_Set_Wgt
(
    FlipsCtx  *FC,
    VSPINT    var,
    VSPINT    wgt
)
{
    SET_FUNCTION_COUNTER ;
    TMP_ASSERT(var == 0 || var == 1) ;
    FC->wgts[var] = wgt ;
}

// ASSUMED: The flip is not empty
// Return 1 for upflip, -1 for downflip
VSPINT FlipsCtx_GetFlipDir
(
    FlipsCtx    *FC,
    Flip      *flip
)
{
    return 1 - 2 * (FC->vars[flip->var])[flip->ind] ;
}

// The flip is acceptable if one of the following conditions is met:
// If this flip is an upflip:
//  1) It is from a better list than the other upflip list, or if it's worse,
//     then there is at least one flip in that list that's infeasible
// 2) If the flip is not from the best linked list, then in all optimal linked
void FlipsCtx_CheckFlip
(
    FlipsCtx    *FC,
    Flip      *flip
)
{

}
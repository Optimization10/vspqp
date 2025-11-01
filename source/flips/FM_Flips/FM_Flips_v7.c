#include "shared_api.h"
#include "FM_Flips/FM_Flips_v7.h"


/*
 * 
 * Lists[v][zi] : Contains a pointer to a Doubly_LLs_Sized struct
 *                corresponding to variable v (0 or 1) with value zi (0 or 1).
 *                For example, Lists[0][0] contains doubly linked lists
 *                associated with x_i = 0
 * 
 * Ts[v][zi] : Contains the optimal Axi or Ayi value corresponding to variable
 *             v (0 or 1) and value zi (0 or 1). For example, 
 *             Ts[0][0] = min{Ay_i : x_i = 0, x_i unlocked}. Positive infinity
 *             (inf of the empty set) is recorded as n+1, and negative infinity
 *             (sup of the empty set) is recorded as -1
 *
 * X, Y : These are states which contain pointers to x, y, Ax, Ay, and weights.
 *        At the beginning, these contain pointers to the original point, but
 *        by the end, they may contain different pointers and values for
 *        weights. However, at the end, the pointers pointers correspond to the
 *        best point 
 */


void FM_Flips_v7
(
    VSP_state              *X,
    VSP_state              *Y,
    VSP_problem      *problem,
    VSPINT       maxNumPasses
)
{
    /*------------------------------------------------------------------------*/
    // Debugging checks for input
    /*------------------------------------------------------------------------*/
    static long long unsigned count = 0 ;
    count++ ;
    #ifdef VSP_DEBUG
        // Check input data
        VSP_check_states(X, Y, problem) ;
        if(maxNumPasses <= 0 ) ERR("Incorrect maxNumPasses") ;

        // Check that x and y are binary
        if(X->point->ind_frac != EMPTY || Y->point->ind_frac != EMPTY)
        {
            ERR("Non-binary points in FM flips") ;
        }
        
        // Record initial objective
        VSPINT objInit ;
        Compute_ObjValSimple
        (
            X->point->int_part, X->Az->int_part, X->weight,
            Y->point->int_part, Y->Az->int_part, Y->weight,
            problem, &objInit
        ) ;
    #endif

    /*------------------------------------------------------------------------*/
    // Encode data in 0-1 format, and as states
    /*------------------------------------------------------------------------*/
    VSPINT *locked[2],
    *vars [2] = { X->point->int_part  ,   Y->point->int_part  },
    *Avars[2] = { X->Az->int_part     ,   Y->Az->int_part     },
    *wgts [2] = { &(X->weight)        ,   &(Y->weight)        }  ;

    VSP_state *retStates[2] = {X, Y} ;


    /*------------------------------------------------------------------------*/
    // Initialize data
    /*------------------------------------------------------------------------*/
    VSPINT k, v, i, indFlip, varFlip, gainFlip,
    nflips = 0,
    n      = problem->Graph->n ;
    // Locked arrays
    for(v = 0 ; v < 2 ; v++) locked[v] = calloc(n, sizeof(VSPINT)) ;

    // Linked Lists
    Doubly_LLs_Sized *Lists[2][2] ;
    for(v = 0 ; v < 2 ; v++)
    {
        for(i = 0 ; i < 2 ; i++)
        {
            Lists[v][i] = Doubly_LLs_Sized_Init(n) ;
        }
    }

    // Ts
    VSPINT Ts[2][2] ;
    Initialize_Lists_And_Ts(vars, Avars, locked, Lists, Ts, problem, TRUE);

    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

    // Memory
    Flip_Record *memory = malloc( 2 * maxNumPasses * n * sizeof(*memory)) ;
    VSPINT wgtImb,
    accumGain     = 0,
    bestAccumGain = 0,
    bestPoint     = 0, // 0 corresponds to initial point, after 0 flips.        
    bestWgtImb    = IMBALANCE(*wgts[0], *wgts[1] ) ;

    // Checkpoints
    VSPINT m = n / 2 ; // Number of flips between checkpoints  

    // H is the largest integer which satisfies Hm <= nflips, where nflips is
    // the total number of flips done. Since the last point will also be a
    // checkpoint, which is not necesarily a multiple of m, the number
    // of checkpoints it genearally H + 1, unless nflips is a multiple of m,
    // in which case the number of checkpoints is H 
    VSPINT H = 0 ;

    VSPINT maxNumCPs = ceil(( ( 2.0*maxNumPasses*n) / m ) + 1) ;
    VSP_state **CPs[2] =
    {
        calloc( maxNumCPs, sizeof(VSP_state*) ),
        calloc( maxNumCPs, sizeof(VSP_state*) )
    } ;

    for(v = 0 ; v < 2 ; v++)
    {
        Save_State_Copy(vars[v], Avars[v], *wgts[v], n, &(CPs[v][0]) ) ;
    }

    /*------------------------------------------------------------------------*/
    // Do Flips
    /*------------------------------------------------------------------------*/
    for(VSPINT numPasses = 0 ; numPasses < maxNumPasses ; numPasses++)
    {
        /*---- Reset and Initialize LLs, Ts, and Locked Arrays ---------------*/
        if(numPasses > 0)
        {
            for(v = 0 ; v < 2 ; v++)
            {
                Locked_Arr_Reset(locked[v], n) ;
                for(i = 0 ; i < 2 ; i++) Doubly_LLs_Sized_Reset(Lists[v][i]) ;
            }
            Initialize_Lists_And_Ts
            (
                vars, Avars, locked, Lists, Ts, problem, FALSE
            ) ;
        }
        IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

        /*---- Flips ---------------------------------------------------------*/
        do
        {
            Select_Flip
            (
                vars     , Avars    , wgts      , locked , Ts , Lists , problem,
                &varFlip , &indFlip , &gainFlip
            ) ;
            IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

            if( varFlip != NO_VAR )
            {
                // Early exit on the last pass
                if(numPasses == maxNumPasses - 1 && gainFlip <= 0 ) break ;

                Do_Flip
                (
                    vars    , Avars   , wgts    , locked, Ts, Lists,
                    varFlip , indFlip , problem
                ) ;

                IFDBG(Check_Data(vars, Avars, wgts, locked, Ts, Lists,problem));
                Update_Local_Data_With_Flip
                (
                    vars, Avars, wgts, m, &nflips, &H, &accumGain, &bestPoint,
                    &bestAccumGain, &bestWgtImb, indFlip, varFlip, gainFlip,
                    CPs, memory, problem
                ) ;
                IFDBG(Check_Data(vars, Avars, wgts, locked, Ts, Lists,problem));
            }
            else break ; // Exit pass if no feasible, unlocked flip
        } while( 1 ) ;
    }

    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

    // If the last point visited was not saved as a checkpoint, do it now.
    if(nflips % m != 0)
    {
        CPs[X_VAR][H+1] = X ;
        CPs[Y_VAR][H+1] = Y ;
    }

    /*------------------------------------------------------------------------*/
    // Rollback
    /*------------------------------------------------------------------------*/
    Rollback(CPs, retStates, m, H, bestPoint, nflips, memory, problem) ;

    /*------------------------------------------------------------------------*/
    // Free memory
    /*------------------------------------------------------------------------*/
    // Free LLs
    for(v = 0 ; v < 2 ; v++)
    {
        for(i = 0 ; i < 2 ; i++) Doubly_LLs_Sized_Destroy(Lists[v][i]) ;
    }

    for(v = 0 ; v < 2 ; v++) free(locked[v]) ;

    // Free memory array
    free(memory) ;

    // Free checkpoints
    for(i = 0 ; i <= H  ; i++)
    {
        for(v = 0 ; v < 2 ; v++) VSP_destruct_state(CPs[v][i]) ;
    }
    for(v = 0 ; v < 2 ; v++) free(CPs[v]) ;

    /*------------------------------------------------------------------------*/
    // Debugging checks for output
    /*------------------------------------------------------------------------*/
    #ifdef VSP_DEBUG
        // Check output data
        VSP_check_states(X, Y, problem) ;

        // Check that x and y are binary
        if(X->point->ind_frac != EMPTY || Y->point->ind_frac != EMPTY)
        {
            ERR("Non-binary points in FM flips") ;
        }
        
        // Record final objective
        VSPINT objFin ;
        Compute_ObjValSimple
        (
            X->point->int_part, X->Az->int_part, X->weight,
            Y->point->int_part, Y->Az->int_part, Y->weight,
            problem, &objFin
        ) ;

        // Check that the objective has improved by bestAccumGain and it has
        // bestImb
        if(objInit + bestAccumGain != objFin)
        {
            ERR("Objective did not increase as promised") ;
        }

        VSPINT imb = IMBALANCE(X->weight, Y->weight) ;
        if(imb != bestWgtImb)
        {
            ERR("Imbalance not as promised") ;
        }

        if(objInit > objFin)
        {
            ERR("Objective got worse in FM Flips") ;
        }
    #endif

    return ;
}


static void Update_Local_Data_With_Flip
(
    VSPINT            *vars[2],
    VSPINT           *Avars[2],
    VSPINT            *wgts[2],
    VSPINT                   m,
    VSPINT             *nflips,
    VSPINT                  *H,
    VSPINT          *accumGain,
    VSPINT          *bestPoint,
    VSPINT      *bestAccumGain,
    VSPINT         *bestWgtImb,
    VSPINT             indFlip,
    VSPINT             varFlip,
    VSPINT            gainFlip,
    VSP_state         **CPs[2],
    Flip_Record        *memory,
    VSP_problem       *problem
)
{
    VSPINT wgtImb, v,
    n = problem->Graph->n ;

    // Save flip in memory
    memory[(*nflips)++] = (Flip_Record) {.var = varFlip, .ind = indFlip} ;

    // Update accumulated gain and best point
    *accumGain += gainFlip ;
    wgtImb     = IMBALANCE(*wgts[X_VAR], *wgts[Y_VAR]) ;
    if(IS_BETTER_POINT(*accumGain,*bestAccumGain,wgtImb,*bestWgtImb))
    {
        *bestPoint     = *nflips ;
        *bestAccumGain = *accumGain ;
        *bestWgtImb    = wgtImb ;
    }

    // Save checkpoint
    if(*nflips % m == 0)
    {
        (*H)++ ;
        for(v = 0 ; v < 2 ; v++)
        {
            Save_State_Copy(vars[v], Avars[v], *wgts[v], n, &(CPs[v][*H]) ) ;
        }
    }
}


static void Select_Flip
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSP_problem          *problem,
    VSPINT            *retFlipVar,
    VSPINT            *retIndFlip,
    VSPINT           *retFlipGain
)
{
    SET_FUNCTION_COUNTER ;
    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

    VSPINT  i, v, zi, k, ci, T, feasGains[4],
    flipInds [4] = {EMPTY, EMPTY, EMPTY, EMPTY},
    numFlips [4] = {0},
    gamma = problem->Graph->max_vtx_cost,
    n     = problem->Graph->n,
    u     = problem->u,
    *c    = problem->Graph->c ;


    /*------------------------------------------------------------------------*/
    // Collect the following data
    // flipInds[4]  - in entry k, record the index of the best flip of type k
    //                  if such a flip exists and is feasible. Otherwise,
    //                  leave the initial value of -1 unchanged
    // feasGains[4] - in entry k, record the feasible gain (if any)
    //                  of the flip encoded as k
    // numFlips[4]  - in entry k, if flip k exists and is feasible, store 
    //                the size of the optimal heap
    // maxFeasGain  - maximal feasible positive gain (weakly positive, >= 0)
    /*------------------------------------------------------------------------*/
    const VSPINT noGainFlag = -gamma * n - gamma  ;
    VSPINT maxFeasGain = noGainFlag ;

    for(v = 0 ; v < 2 ; v++)
    {
        for(zi = 0 ; zi < 2 ; zi++)
        {
            T = Ts[v][zi] ;
            if( IS_VALID_T(zi, T) )
            {
                Doubly_LLs_Sized_See_Entry(Lists[v][zi], T, &i) ;
                VSP_ASSERT(i != EMPTY) ;

                // If flip is feasible, compute gain and store it in feasGains
                ci = c[i] ;
                if( IS_FEASIBLE_FLIP(zi, *wgts[v], ci, u ) )
                {
                    k = GET_ENCODING(v, zi) ;
                    flipInds [k] = i ;
                    numFlips [k] = Lists[v][zi]->sizes[T] ;
                    feasGains[k] = GAIN_FLIP( DIR_FLIP(zi), ci, gamma, T ) ;
                    if( feasGains[k] > maxFeasGain )
                    {
                        maxFeasGain = feasGains[k] ;
                    } 
                }
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Quit early if there is no feasible gain
    /*------------------------------------------------------------------------*/
    if(maxFeasGain == noGainFlag)
    {
        *retFlipVar = NO_VAR ;
        *retIndFlip = EMPTY ;
        *retFlipGain = maxFeasGain ;
        IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
        return ;
    }

    /*------------------------------------------------------------------------*/
    // Record all flips with maximal gain in 
    // maxGainFlips[4]     - a stack containing 0-3 which encode flips
    // m                   - size of maxGainFlips stack
    // numsMaxGainFlips[4] - a stack of the number of indices corr. to each flip
    // accum[5]            - running total of numsMaxGainFlips
    /*------------------------------------------------------------------------*/
    VSPINT maxGainFlips[4]     = {0} ;
    VSPINT m                   =  0  ;
    VSPINT numsMaxGainFlips[4] = {0} ;
    VSPINT acumm[5]            = {0} ;

    for(k = 0 ; k < 4 ; k++)
    {
        if(flipInds[k] != EMPTY && feasGains[k] == maxFeasGain)
        {
            maxGainFlips    [m] = k ;
            numsMaxGainFlips[m] = numFlips[k] ;
            m++ ; 
            acumm[m]            = acumm[m-1] + numFlips[k] ;
        }
    }
    VSP_ASSERT(m >= 1) ;


    /*------------------------------------------------------------------------*/
    // Choose a random flip with maximal gain
    // flipChoice - integer in [0, 3] representing the flip
    /*------------------------------------------------------------------------*/
    VSPINT randChoice  = rand() % acumm[m] ;
    VSPINT flipChoice ;


    for(VSPINT p = 1 ; p <= m ; p++)
    {
        if(randChoice < acumm[p])
        {
            flipChoice = maxGainFlips[p-1] ;
            break ;
        }
    }
    *retFlipVar  = GET_VAR(flipChoice) ;
    *retIndFlip  = flipInds[flipChoice] ;
    *retFlipGain = maxFeasGain ;

    VSPINT dir = GET_DIR(flipChoice) ;

    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
    return ;
}




// Let z be the variable of the flip, and w the other variable
// ASSUMED: The flip is happening with an Awi or Azi with optimal T
static void Do_Flip
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSPINT                    var, // variable of the flip
    VSPINT                      i, // index of the flip
    VSP_problem          *problem
)
{
    SET_FUNCTION_COUNTER ;
    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;

    VSPLONG p ;
    VSPINT j, T, wj, Azj,
    otherVar = OTHER_VAR(var) ;

    VSPINT
    n   = problem->Graph->n,
    *z  = vars [var],
    *Az = Avars[var],
    *w  = vars [otherVar],
    *Aw = Avars[otherVar],
    zi  = vars [var][i] ;

    VSPINT
    dir = DIR_FLIP(zi),
    Awi = Ts[var][zi],
    *c  = problem->Graph->c,
    *Ai = problem->Graph->Ai ;
    VSPLONG *Ap = problem->Graph->Ap ;

    Doubly_LLs_Sized *lists ;

    /*------------------------------------------------------------------------*/
    // Remove i from its list
    // Udate the T corresponding to zi
    // Update weights
    // Update zi to di + dir and lock zi
    /*------------------------------------------------------------------------*/
    Doubly_LLs_Sized_Delete(Lists[var][zi], Awi, i) ;

    // If i was the only entry in its list, then the T needs to be updated
    T = Awi ;
    if( zi ) // zi = 1
    {
        // We are interested in max, so search downwards
        while(T > -1 && Lists[var][zi]->sizes[T] == 0 ) T-- ;
    }
    else
    {
        // We are interested in min, so search upwards
        while(T < n + 1 && Lists[var][zi]->sizes[T] == 0 ) T++ ;
    }

    // Update T
    Ts[var][zi] = T ;

    // Update zi, weight, and locked status
    z[i]       += dir ;
    *wgts[var] += dir * c[i] ;
    locked[var][i]++ ;


    /*------------------------------------------------------------------------*/
    // Loop through j in N(i) cup (i) where wj is unlocked, and shift them one
    // list up or down
    /*------------------------------------------------------------------------*/
    VSPLONG start = Ap[i] ;
    VSPLONG end   = Ap[i+1] ;
    p = start ;
    j = i ;

    bool iItselfFinished = false ;
    do
    {
        if(iItselfFinished)
        {
            j = Ai[p] ;
            p++ ;
        }
        iItselfFinished = true ;

        // if wj is locked, do nothing
        if( !locked[otherVar][j] )
        {
            // Remove j from its list
            wj    = w[j] ;
            Azj   = Az[j] ;
            lists = Lists[otherVar][wj] ;

            Doubly_LLs_Sized_Delete(lists, Azj, j) ;
            Doubly_LLs_Sized_Insert(lists, Azj + dir, j) ;

            // If wj = 0 and dir = 1, or if wj = 1 and dir = -1
            if(dir + 2 * wj == 1)
            {
                T = Ts[otherVar][wj] ;
                if(lists->sizes[T] == 0 )
                {
                    Ts[otherVar][wj] += dir ;
                }
            }
            else // If wj = 0 and dir = -1, or if wj = 1 and dir = 1
            {
                if( Azj == Ts[otherVar][wj])
                {
                    Ts[otherVar][wj] += dir ;
                }
            }     
        }
        Az[j] += dir ;
    } while ( p != end ) ;

    IFDBG( Check_Data(vars, Avars, wgts, locked, Ts, Lists, problem) ) ;
    return ;
}


static void Rollback
(
    VSP_state            **CPs[2],
    VSP_state       *retStates[2],
    VSPINT                      m,
    VSPINT                      H,
    VSPINT              bestPoint,
    VSPINT                 nflips,
    Flip_Record           *memory,
    VSP_problem          *problem
)
{
    
    VSPINT i, v, p, nflipsFromLo, nflipsFromHi, whichCP, nflipsToDo,
    start, // an index of mem containing the first flips which should be done
    dir,   // dir - +1 if we go from lower cp up, and -1 if from higher cp down
    *c = problem->Graph->c,
    *Ai = problem->Graph->Ai ;
    VSPLONG *Ap = problem->Graph->Ap ;

    VSPINT h ; // h has the property that hm <= bestPoint < (h+1)m

    VSP_state *CPChosen[2] ;
    Flip_Record flip ;

    // Determine the closest checkpoint to bestPoint
    if(bestPoint < H * m)
    {
        h = (bestPoint - bestPoint % m) / m ;
        nflipsFromLo = bestPoint - h * m ;
        nflipsFromHi = (h + 1) * m  - bestPoint ;
        whichCP      = (nflipsFromLo <= nflipsFromHi) ? h : h + 1 ;
        start        = (nflipsFromLo <= nflipsFromHi) ? h * m : (h + 1) * m - 1;
    }
    else
    {
        nflipsFromLo = bestPoint - H * m ;
        nflipsFromHi = nflips - bestPoint ;
        // The <= catches the case when bestPoint = nflips and nflips % m = 0
        // In this case, there is no (H+1)th checkpoint, so we MUST "roll" from
        // the lower checkpoint
        whichCP      = (nflipsFromLo <= nflipsFromHi) ? H : H + 1 ;
        start        = (nflipsFromLo <= nflipsFromHi) ? H * m : nflips - 1 ;
    }

    for(v = 0 ; v < 2 ; v++)
    {
        CPChosen[v] = CPs[v][whichCP] ;
    }
    nflipsToDo   = MIN(nflipsFromLo, nflipsFromHi) ;
    dir          = (nflipsFromLo <= nflipsFromHi) ? 1 : -1 ;


    for(p = 0 ; p < nflipsToDo ; p++)
    {
        flip = memory[start + dir * p] ;
        v = flip.var ;
        i = flip.ind ;
        Do_Flip_Basic
        (
            CPChosen[v]->point->int_part  ,  CPChosen[v]->Az->int_part  ,   i,
            &(CPChosen[v]->weight)        ,  c, Ai                ,   Ap
        ) ;
    }

    // Swap contents of cp and return states
    for(v = 0 ; v < 2 ; v++)
    {
        Swap_Mem(CPChosen[v], retStates[v], sizeof(VSP_state) ) ;
    }

    return ;
}

static void Locked_Arr_Reset(VSPINT *locked, VSPINT n)
{
    SET_FUNCTION_COUNTER ;
    VSP_ASSERT(locked != NULL) ;
    for(VSPINT i = 0 ; i < n ; i++) locked[i] = 0 ;
    return ;
}


// Initialize sized doubly linked lists and Ts. Ignore locked vertices. If
// randomize is true, then build the lists in a randomized manner by going
// through the vertices in a random order. Otherwise, go through the vertices
// in numerical order
// locked - If null, then ignored. One or more of its entries can also be NULL,
//          in which case that particular array is ignored
void Initialize_Lists_And_Ts
(
    VSPINT                 *vars[2],
    VSPINT                *Avars[2],
    VSPINT               *locked[2],
    Doubly_LLs_Sized   *lists[2][2],
    VSPINT                 Ts[2][2],
    VSP_problem            *problem,
    VSPINT                randomize  // 0 or 1
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i, v, j,*perm, zi, Awi, encoding, dir,
    n = problem->Graph->n ;
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
        Ts[v][0] = n + 1 ;
        Ts[v][1] = -1 ;
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
            Doubly_LLs_Sized_Insert(lists[v][zi], Awi, i);
            if(dir * Awi < dir * Ts[v][zi] ) Ts[v][zi] = Awi ;
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



static void Check_Data
(
    VSPINT               *vars[2],
    VSPINT              *Avars[2],
    VSPINT               *wgts[2],
    VSPINT             *locked[2],
    VSPINT               Ts[2][2],
    Doubly_LLs_Sized *Lists[2][2],
    VSP_problem          *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT v, zi, T, Tdbg, k,
    n = problem->Graph->n ;
    /*------------------------------------------------------------------------*/
    // Check vars, Avars, weights, and feasibility by placing into states
    /*------------------------------------------------------------------------*/
    VSP_point
    ptX = {.int_part = vars[X_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n},
    ptY = {.int_part = vars[Y_VAR], .ind_frac = EMPTY, .val_frac = 0, .n = n} ;

    VSP_Az
    structAx = {.int_part   = Avars[X_VAR] , .ind_frac_z = EMPTY , .n = n   ,
                .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   },
    structAy = {.int_part   = Avars[Y_VAR] , .ind_frac_z = EMPTY , .n = n   ,
                .val_frac_z = 0            , .col_frac_z = NULL  , .one = 1   };

    VSP_state
    X = {.point = &ptX, .Az = &structAx, .weight = *wgts[X_VAR]},
    Y = {.point = &ptY, .Az = &structAy, .weight = *wgts[Y_VAR]} ;
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
        vars[X_VAR], Avars[X_VAR], *wgts[X_VAR], locked[X_VAR],
        vars[Y_VAR], Avars[Y_VAR], *wgts[Y_VAR], locked[Y_VAR],
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
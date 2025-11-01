#include "FM_Flips_Internal.h"
#include "ds_api.h"
#include "limits.h"
#include "stdlib.h"
#include "time.h"

extern double  G_ROLLBACK_TIME ;

#define TEMP

#define IS_BETTER_POINT(accumGain, bestAccumGain, wgtImb, bestWgtImb) \
  ((accumGain) > (bestAccumGain) || \
  ((accumGain) == (bestAccumGain) && (wgtImb) <= (bestWgtImb)))



// 1) This function picks a random best flip (but does not apply it)
Flip Get_Best_Flip_Random
(
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;

// 2) This function applies (executes) the flip returned by Get_Best_Flip_Random
void Do_Given_Flip
(
    Flip            flip,
    FM_Flips_Struct *zstruct,
    FM_Flips_Struct *wstruct,
    VSP_problem     *problem
) ;


void FM_Flips_API_v2
(
    VSP_state         *Z,
    VSP_state         *W,
    VSP_problem *problem,
    VSP_options *options
)
{
    TIMER_START ;
    SET_FUNCTION_COUNTER ;

    #ifdef VSP_DEBUG
        VSP_check_states(Z, W, problem) ;
        VSPINT nFlipsLastPass = 0 ;   
    #endif

    struct timespec startT, endT ;
    double elapsed ;

    // Local vars
    VSPINT n       = problem->Graph->n,
           gamma   = problem->Graph->max_vtx_cost,
           wgtzOrg = Z->weight,
           wgtwOrg = W->weight,
          *Ai      = problem->Graph->Ai,
          *Ap      = problem->Graph->Ap,
          *c       = problem->Graph->c,
          *z       = Z->point->int_part,
          *Az      = Z->Az->int_part,
          *w       = W->point->int_part,
          *Aw      = W->Az->int_part,
          i, p ;
    Flip *memory = malloc (8 * n * sizeof(Flip)) ;

    #ifdef VSP_DEBUG
        VSPFLOAT objin = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
    #endif

    // Initialize the structs
    FM_Flips_Struct *Z_struct = Construct_FM_Flips_Struct
                                (
                                    z, Az, Aw, c, n, wgtzOrg, gamma, Ai, Ap, problem->u
                                ) ;
    Z_struct->var = 'z' ;

    FM_Flips_Struct *W_struct = Construct_FM_Flips_Struct
                                (
                                    w, Aw, Az, c, n, wgtwOrg, gamma, Ai, Ap, problem->u
                                ) ;
    W_struct->var = 'w' ;


    VSPINT nflips = 0 ;
    Flip flip ;
    VSPINT accumGain       = 0,
           bestAccumGain   = 0,
           bestPoint       = 0, // 0 corresponds to initial point, after 0 flips.        
           bestWgtImb      = (VSPINT) abs(Z_struct->weight - W_struct->weight),
           maxNumPasses    = 4,
          *numFlipsPerPass = calloc(maxNumPasses + 1, sizeof(VSPINT)),
           passBestPoint   = EMPTY,
           wgtImb ;

    // printf("\n\n\nStarting New pass") ;
    /*------------------------------------------------------------------------*/
    // Variables for checkpoints
    /*------------------------------------------------------------------------*/
    // Background: The flips are indexed starting at 1.
    // Each time we do a flip, we go from one point to the next.
    // The points we visit are indexed starting at 0. Hence, the original point
    // is point 0, the point after the first flip is point 1, etc.
    // The original point is the first checkpoint, and after every m flips, we
    // save a new checkpoint. The checkpoints are also indexed starting at 0.
    // The final point reached is also considered a checkpoint, even though
    // there may be less than m flips between it and the last checkpoint.
    // For this last point, we do not save a copy of the point, but rather
    // use input variables themselves. However, this rule is violated if
    // the total number of flips is a multiple of m. The reason is that the
    // decision to terminate flips on the last pass of the loop is only decided
    // in the iteration of the loop that comes after the iteration in which we
    // save the checkpoint

    // A positive integer. Number of flips between checkpoints
    // VSPINT m = n / 4 ;
    VSPINT m = n / 2 ;
    // VSPINT m = 2 * n / 3 ;
    // VSPINT m = n ;
    // VSPINT m = (VSPINT) (1.5 * n) ;
    // VSPINT m = 2 * n ;
    // VSPINT m = 2 * n ;
    // VSPINT m = 9 * n ;

    // h has the property that hm <= bestPoint < (h+1)m
    VSPINT h ;

    // H is the largest integer which satisfies Hm <= nflips, where nflips is
    // the total number of flips done. Since the last point will also be a
    // checkpoint, which is not necesarily a multiple of m, the number
    // of checkpoints it genearally H + 1, unless nflips is a multiple of m,
    // in which case the number of checkpoints is H 
    VSPINT H = 0 ;

    // Checkpoint p is stored in yCPs
    clock_gettime(CLOCK_MONOTONIC, &startT) ;
    VSP_state **xCPs = calloc( ceil(( (8.0 * n) / m ) + 1), sizeof(VSP_state*) ),
              **yCPs = calloc( ceil(( (8.0 * n) / m ) + 1), sizeof(VSP_state*) );
    Save_State_Copy( z, Az, Z->weight, n, &(xCPs[0]) ) ;
    Save_State_Copy( w, Aw, W->weight, n, &(yCPs[0]) ) ;
    clock_gettime(CLOCK_MONOTONIC, &endT) ;
    elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9;
    G_ROLLBACK_TIME += elapsed ;

    for(VSPINT numPasses = 0 ; numPasses < maxNumPasses ; numPasses++)
    {
        FM_Flips_Struct_Build(Z_struct, problem->Graph->max_vtx_cost) ;
        FM_Flips_Struct_Build(W_struct, problem->Graph->max_vtx_cost) ;
        do
        {
            #ifdef VSP_DEBUG
                FM_Flips_Struct_Check_v2(Z_struct, true) ;
                FM_Flips_Struct_Check_v2(W_struct, true) ;
            #endif

            flip = Get_Best_Flip_Random(Z_struct, W_struct, problem) ;

            if(flip.var != 'n')
            {
                if(numPasses == maxNumPasses - 1 && flip.gain <= 0) break ;
                Do_Given_Flip(flip, Z_struct, W_struct, problem) ;
                memory[nflips++] = flip ;
                accumGain += flip.gain ;
                wgtImb     = abs(Z_struct->weight - W_struct->weight) ;

                // Save checkpoint
                clock_gettime(CLOCK_MONOTONIC, &startT) ;
                if(nflips % m == 0)
                {
                    H++ ;
                    Save_State_Copy( z, Az, Z_struct->weight, n, &(xCPs[H]) ) ;
                    Save_State_Copy( w, Aw, W_struct->weight, n, &(yCPs[H]) ) ;
                    #ifdef VSP_DEBUG
                        VSP_check_states(xCPs[H], yCPs[H], problem) ;
                    #endif
                }
                clock_gettime(CLOCK_MONOTONIC, &endT) ;
                elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9;
                G_ROLLBACK_TIME += elapsed ;

                if(IS_BETTER_POINT(accumGain,bestAccumGain,wgtImb,bestWgtImb))
                {
                    bestPoint     = nflips ;
                    bestAccumGain = accumGain ;
                    passBestPoint = numPasses ;
                    bestWgtImb    = wgtImb ;
                }
            }
            else
            {
                FM_Flips_Struct_Reset(Z_struct) ;
                FM_Flips_Struct_Reset(W_struct) ;
                numFlipsPerPass[numPasses+1] = nflips ;
                break ;
            }
        } while( true ) ;
    }

    // If the last point visited was not saved as a checkpoint, do it now.
    // clock_gettime(CLOCK_MONOTONIC, &startT) ;
    Z->weight = Z_struct->weight ;
    W->weight = W_struct->weight ;
    if(nflips % m != 0)
    {
        xCPs[H+1] = Z ;
        yCPs[H+1] = W ;

    }
    // clock_gettime(CLOCK_MONOTONIC, &endT) ;
    // elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9;
    // G_ROLLBACK_TIME += elapsed ;



    /*------------------------------------------------------------------------*/
    // Roll back
    /*------------------------------------------------------------------------*/
    clock_gettime(CLOCK_MONOTONIC, &startT) ;
    // start - an index of mem containing the first flips which should be done
    VSPINT start ;
    // dir - +1 if we go from lower cp up, and -1 if from higher cp down
    VSPINT dir ;
    VSPINT nflipsFromLo, nflipsFromHi, whichCP, nflipsToDo ;
    VSP_state *xCP, *yCP ;

    // Determine the closest checkpoint to bestPoint
    if(bestPoint < H * m)
    {
        h = (bestPoint - bestPoint % m) / m ;
        nflipsFromLo = bestPoint - h * m ;
        nflipsFromHi = (h + 1) * m  - bestPoint ;
        whichCP      = (nflipsFromLo <= nflipsFromHi) ? h : h + 1 ;
        start        = (nflipsFromLo <= nflipsFromHi) ? h * m : (h + 1) * m - 1 ;
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
    xCP          = xCPs[whichCP] ;
    yCP          = yCPs[whichCP] ;
    nflipsToDo   = MIN(nflipsFromLo, nflipsFromHi) ;
    dir          = (nflipsFromLo <= nflipsFromHi) ? 1 : -1 ;

    #ifdef VSP_DEBUG
        VSP_check_states(xCP, yCP, problem) ;
    #endif
    for(p = 0 ; p < nflipsToDo ; p++)
    {
        flip = memory[start + dir * p] ;
        if(flip.var == 'z')
            Do_Flip_Basic
            (
                xCP->point->int_part, xCP->Az->int_part, flip.ind, &(xCP->weight) ,c, Ai, Ap
            ) ;
        else
        {
            Do_Flip_Basic
            (
                yCP->point->int_part, yCP->Az->int_part, flip.ind, &(yCP->weight) ,c, Ai, Ap
            ) ;
        } 
        #ifdef VSP_DEBUG
            VSP_check_states(xCP, yCP, problem) ;
        #endif
    }


    // Switch contents of cp and return states
    Swap_Mem(xCP, Z, sizeof(VSP_state) ) ;
    Swap_Mem(yCP, W, sizeof(VSP_state) ) ;

    clock_gettime(CLOCK_MONOTONIC, &endT) ;
    elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9;
    G_ROLLBACK_TIME += elapsed ;
    /*------------------------------------------------------------------------*/
    // Old rollback
    /*------------------------------------------------------------------------*/
    // clock_gettime(CLOCK_MONOTONIC, &startT) ;

    // for(VSPINT p = nflips - 1 ; p >= bestPoint ; p--)
    // {
    //     flip     = memory[p] ;
    //     char var = flip.var ;
    //     VSPINT i = flip.ind ;

    //     VSPINT dir ;

    //     #ifdef VSP_DEBUG
    //         VSPFLOAT obj1 = Objective_No_Wgt_Check(Z, W, problem, problem->Graph->max_vtx_cost) ;
    //     #endif

    //     if(var == 'z')
    //     {
    //         dir = 1 - 2 * Z->point->int_part[i] ;
    //         Z->point->int_part[i] += dir ;
    //         Z_struct->weight += (dir * problem->Graph->c[i]) ;
    //         VSP_Az_pm_col_of_A(Z->Az->int_part, i, dir, problem) ;
    //     }
    //     else if(var == 'w')
    //     {
    //         dir = 1 - 2 * W->point->int_part[i] ;
    //         W->point->int_part[i] += dir ;
    //         W_struct->weight += (dir * problem->Graph->c[i]) ;
    //         VSP_Az_pm_col_of_A(W->Az->int_part, i, dir, problem) ;
    //     }
    //     else
    //     {
    //         #ifdef VSP_DEBUG
    //             VSP_error("Flip in non-accepted variable", DBG_INFO) ;
    //         #endif
    //     }

    //     Z->weight = Z_struct->weight ;
    //     W->weight = W_struct->weight ;

    //     #ifdef VSP_DEBUG
    //         VSP_check_states(Z, W, problem) ;
    //         VSPFLOAT obj2 = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;

    //         // Negative flip accumGain because going backwards
    //         VSPINT diffobj = obj2 - obj1 ;
    //         UNUSED_VAR(diffobj) ;
    //         if(obj2 - obj1 != -flip.gain)
    //         {
    //             VSP_error("Gain incorrect", DBG_INFO) ;
    //         }
    //     #endif

    // }
    // Z->weight = Z_struct->weight ;
    // W->weight = W_struct->weight ;

    // clock_gettime(CLOCK_MONOTONIC, &endT);
    // elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9 ;
    // G_ROLLBACK_TIME += elapsed ;

    /*------------------------------------------------------------------------*/
    // Free memory
    /*------------------------------------------------------------------------*/
    free(memory) ;
    Destruct_FM_Flips_Struct(Z_struct) ;
    Destruct_FM_Flips_Struct(W_struct) ;
    free(numFlipsPerPass) ;
    // clock_gettime(CLOCK_MONOTONIC, &startT) ;
    for(VSPINT i = 0 ; i <= H  ; i++)
    {
        VSP_destruct_state(xCPs[i]) ;
        VSP_destruct_state(yCPs[i]) ;
    }
    // if(nflips % m != 0)
    // {
    //     VSP_destruct_state(xCPs[H+1]) ;
    //     VSP_destruct_state(yCPs[H+1]) ;
    // }
    free(xCPs) ;
    free(yCPs) ;
    // clock_gettime(CLOCK_MONOTONIC, &endT) ;
    // elapsed = (endT.tv_sec - startT.tv_sec) + (endT.tv_nsec - startT.tv_nsec) / 1e9;
    // G_ROLLBACK_TIME += elapsed ;

    #ifdef VSP_DEBUG
        VSP_check_states(Z, W, problem) ;
        VSPFLOAT objout = VSP_compute_obj(Z, W, problem, problem->Graph->max_vtx_cost) ;
        VSPFLOAT dobj = objout - objin ;
        if(dobj < 0)
        {
            VSP_error("FM_Flips worsened the objective", DBG_INFO) ;
        }

        if(dobj != bestAccumGain)
        {
            VSP_error("Change in objective not as predicted", DBG_INFO) ; 
        }
    #endif

    TIMER_END ;
    PERFORMANCE_UPDATE_TIME(FM_Flips) ;
}

// numFlips contains non-negative integers. Ex {0, 0, 5, 3}
// This function would return (index) 2 with probability 5/8, and (index) 3
// with probability 3/8
VSPINT SelectBestFlipRandomly
(
    VSPINT numFlips[4] 
)
{
    VSPINT cummNumFlips[4] ;
    VSPINT p, sum = 0 ;
    for(p = 0 ; p < 4 ; p++)
    {
        sum += numFlips[p] ;
        cummNumFlips[p] = sum ;
    }

    VSPINT randInt = rand() % sum ;

    // Find the smallest index in cummNumFlips where the random VSPINT is surpassed
    for(p = 0 ; p < 4 ; p++)
    {
        if(cummNumFlips[p] > randInt)
        {
            break ;
        }
    }
    return p ;
}


Flip Get_Best_Flip_Random(FM_Flips_Struct *zstruct, FM_Flips_Struct *wstruct, VSP_problem *problem)
{
    VSPINT n        = problem->Graph->n,
           gamma    = problem->Graph->max_vtx_cost,
           u        = problem->u,
           rand_int = rand();

    VSPINT wgtz = zstruct->weight, wgtw = wstruct->weight;

    Flip  ret;
    VSPINT   d_obj[4], ind_flips[4] = {EMPTY, EMPTY, EMPTY, EMPTY};
    VSPINT   best_flip = EMPTY;
    VSPINT best_gain = -n*gamma - gamma;

    // z down
    if(zstruct->p_hi != -1)
    {
        d_obj[ZDN] = -(zstruct->c_lo - gamma*zstruct->p_hi);
        ind_flips[ZDN] = zstruct->O_dn->list[rand_int % zstruct->O_dn->size] ;
    }
    // w down
    if(wstruct->p_hi != -1)
    {
        d_obj[WDN] = -(wstruct->c_lo - gamma*wstruct->p_hi);
        ind_flips[WDN] = wstruct->O_dn->list[rand_int % wstruct->O_dn->size];
    }
    // z up
    if(zstruct->p_lo != n+1 && (wgtz + zstruct->c_hi <= u))
    {
        d_obj[ZUP] = zstruct->c_hi - gamma*zstruct->p_lo;
        ind_flips[ZUP] = zstruct->O_up->list[rand_int % zstruct->O_up->size];
    }
    // w up
    if(wstruct->p_lo != n+1 && (wgtw + wstruct->c_hi <= u))
    {
        d_obj[WUP] = wstruct->c_hi - gamma*wstruct->p_lo;
        ind_flips[WUP] = wstruct->O_up->list[rand_int % wstruct->O_up->size];
    }

    for(VSPINT p = 0; p < 4; p++)
    {
        if(ind_flips[p] != EMPTY && d_obj[p] > best_gain)
        {
            best_gain = d_obj[p];
            best_flip = p;
        }
    }
    if(best_flip == EMPTY)
    {
        ret.ind = EMPTY; ret.var = 'n'; ret.gain = 0 ; 
        return ret ;
    }

    // Check ties
    VSPINT cnt[4] = {0,0,0,0};
    VSPINT size[4] = { zstruct->O_dn->size, wstruct->O_dn->size,
                       zstruct->O_up->size, wstruct->O_up->size };

    for(VSPINT p = 0; p < 4; p++)
        if(ind_flips[p] != EMPTY && d_obj[p] == best_gain) cnt[p] = size[p] ;

    best_flip = SelectBestFlipRandomly(cnt);
    ret.ind   = ind_flips[best_flip];
    ret.gain  = best_gain;
    ret.var   = (best_flip == ZDN || best_flip == ZUP) ? 'z' : 'w';
    return ret;
}

void Do_Given_Flip(Flip flip, FM_Flips_Struct *zstruct, FM_Flips_Struct *wstruct, VSP_problem *problem)
{
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(zstruct, true) ;
        FM_Flips_Struct_Check_v2(wstruct, true) ;
    #endif
    if(flip.ind == EMPTY) return;
    VSPINT i = flip.ind;
    VSPINT *Ai = problem->Graph->Ai, *Ap = problem->Graph->Ap;
    VSPINT *Az = wstruct->Aw, *Aw = zstruct->Aw;
    VSPINT *z  = zstruct->z,  *w  = wstruct->z;

    if(flip.var == 'z')
    {
        if(z[i])
        {
            FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v2(zstruct, i, zstruct->O_dn->lookup[i]);
            for(VSPLONG p = Ap[i]; p < Ap[i+1]; p++)
            {
                VSPINT j = Ai[p];
                if(w[j]) FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2(wstruct, j, Az[j]);
                else     FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2(wstruct, j, Az[j]);
            }
            if(w[i]) FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2(wstruct, i, Az[i]);
            else     FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2(wstruct, i, Az[i]);
        } else
        {
            FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v2(zstruct, i, zstruct->O_up->lookup[i]);
            for(VSPLONG p = Ap[i]; p < Ap[i+1]; p++)
            {
                VSPINT j = Ai[p];
                if(w[j]) FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2(wstruct, j, Az[j]);
                else     FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2(wstruct, j, Az[j]);
            }
            if(w[i]) FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2(wstruct, i, Az[i]);
            else     FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2(wstruct, i, Az[i]);
        }
    } else {
        if(w[i]) {
            FM_Flips_Struct_Update_zi_1_To_0_Where_i_Is_In_O_dn_v2(wstruct, i, wstruct->O_dn->lookup[i]);
            for(VSPLONG p = Ap[i]; p < Ap[i+1]; p++) {
                VSPINT j = Ai[p];
                if(z[j]) FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2(zstruct, j, Aw[j]);
                else     FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2(zstruct, j, Aw[j]);
            }
            if(z[i]) FM_Flips_Struct_Update_wj_1_And_Azj_Decreases_By_1_v2(zstruct, i, Aw[i]);
            else     FM_Flips_Struct_Update_wj_0_And_Azj_Decreases_By_1_v2(zstruct, i, Aw[i]);
        } else {
            FM_Flips_Struct_Update_zi_0_To_1_Where_i_Is_In_O_up_v2(wstruct, i, wstruct->O_up->lookup[i]);
            for(VSPLONG p = Ap[i]; p < Ap[i+1]; p++) {
                VSPINT j = Ai[p];
                if(z[j]) FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2(zstruct, j, Aw[j]);
                else     FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2(zstruct, j, Aw[j]);
            }
            if(z[i]) FM_Flips_Struct_Update_wj_1_And_Azj_Increases_By_1_v2(zstruct, i, Aw[i]);
            else     FM_Flips_Struct_Update_wj_0_And_Azj_Increases_By_1_v2(zstruct, i, Aw[i]);
        }
    }
    #ifdef VSP_DEBUG
        FM_Flips_Struct_Check_v2(zstruct, true) ;
        FM_Flips_Struct_Check_v2(wstruct, true) ;
    #endif
}
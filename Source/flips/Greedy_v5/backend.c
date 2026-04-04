#include "core.h"
#include "flips_private_v5.h"

void z0Aw1_List_Check_v5
(
    Direct_Access_List *List,
    VSPINT                *c,
    VSPINT                *z,
    VSPINT               *Aw,
    VSPINT             gamma,
    VSPINT                 n
)
{
    SET_FUNCTION_COUNTER ;
    WARN_IF_NOT_DEBUG_MODE ;

    #ifdef DEBUG
    Direct_Access_List_Check(List) ;
    VSPINT i, p ;
    VSPINT *list = List->list ;
    VSPINT size  = List->size ;
    
    // Get a witness array of all i in the list, and make sure all i satisfy
    // zi = 0, Awi = 1, ci = gamma
    VSPINT *witnessList = calloc(n, sizeof(VSPINT)) ;
    for(p = 0 ; p < size ; p++)
    {
        i = list[p] ;
        witnessList[i]++ ;
        if(z[i] != 0)
        {
            ERR("i with zi = 1 in z0Aw1 list") ;
        }

        if(Aw[i] != 1)
        {
            ERR("i with Awi != 1 in z0Aw1 list") ;
        }

        if(c[i] != gamma)
        {
            ERR("i with  ci != gamma in z0Aw1 list") ;
        }
    }

    VSPINT *witnessAccurate = calloc(n, sizeof(VSPINT)) ;
    for(i = 0 ; i < n ; i++)
    {
        if(z[i] == 0 && Aw[i] == 1 && c[i] == gamma)
        {
            witnessAccurate[i]++ ;
        }
    }

    for(i = 0 ; i < n ; i++)
    {
        if(witnessAccurate[i] != witnessList[i])
        {
            ERR("Mismatch between witnesses") ;
        }

        if(witnessList[i] != 0 && witnessList[i] != 1)
        {
            ERR("An i is in a list more than once") ;
        }

        if(witnessAccurate[i] != 0 && witnessAccurate[i] != 1)
        {
            ERR("Very strange error") ;
        }
    }


    free(witnessList) ;
    free(witnessAccurate) ;
    #endif
    return ;
}

// Assume that Upflip_v5 happens in z
void Upflip_v5
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtz,
    VSPINT                       *wgtw,
    IdxGroupsByC         *z0Aw0_groups,
    Direct_Access_List     *z0Aw1_list,
    IdxGroupsByC         *z1Aw1_groups,
    IdxGroupsByC         *z1Aw2_groups,

    IdxGroupsByC         *w0Az0_groups,
    Direct_Access_List     *w0Az1_list,
    IdxGroupsByC         *w1Az1_groups,
    IdxGroupsByC         *w1Az2_groups,

    VSPINT                     indFlip,
    Problem               *problem
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT n = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c = problem->graph->c ;
    VSPINT *Ai   = problem->graph->Ai ;
    VSPLONG *Ap   = problem->graph->Ap ;
    VSPINT i, k, Azk ;
    VSPLONG start, end, p ;

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif



    z[indFlip] = 1 ;
    *wgtz += c[indFlip] ;

    /*------------------------------------------------------------------------*/
    // Do all changes pertaining to indFlip leaving a struct
    /*------------------------------------------------------------------------*/
    VSPINT AwindFlip = Aw[indFlip] ;
    if( AwindFlip == 0 )
    {
        IdxGroupsByC_Delete(z0Aw0_groups, indFlip) ;
    }
    else if( AwindFlip == 1 && c[indFlip] == gamma)
    {
        Direct_Access_List_Delete_By_Value(z0Aw1_list, indFlip) ;
    }

    /*------------------------------------------------------------------------*/
    // Do all changes pertaining to indFlip entering a struct
    /*------------------------------------------------------------------------*/
    if(AwindFlip == 1)
    {
        IdxGroupsByC_Insert(z1Aw1_groups, indFlip) ;
    }
    else if(AwindFlip == 2)
    {
        IdxGroupsByC_Insert(z1Aw2_groups, indFlip) ;
    }


    /*------------------------------------------------------------------------*/
    // Do all changes pertaining Awk increasing in N(indFlip) \cup {indFlip}
    /*------------------------------------------------------------------------*/
    start = Ap[indFlip] ;
    end   = Ap[indFlip + 1] ;
    for(p = start ; p < end ; p++)
    {
        k   = Ai[p] ;
        Azk = Az[k]++ ;
        if(w[k]) // wk = 1
        {
            if(Azk == 0)
            {
                IdxGroupsByC_Insert(w1Az1_groups, k) ;
            }
            else if(Azk == 1)
            {
                IdxGroupsByC_Delete(w1Az1_groups, k) ;
                IdxGroupsByC_Insert(w1Az2_groups, k) ;
            }
            else if(Azk == 2)
            {
                IdxGroupsByC_Delete(w1Az2_groups, k) ;
            }
        }
        else // wk = 0
        {
            if(Azk == 0)
            {
                IdxGroupsByC_Delete(w0Az0_groups, k) ;
                if(c[k] == gamma)
                {
                    Direct_Access_List_Insert(w0Az1_list, k) ;
                }
            }
            else if(Azk == 1)
            {
                if(c[k] == gamma)
                {
                    Direct_Access_List_Delete_By_Value(w0Az1_list, k) ;
                }
            }
        }
    }

    /*------------------------------------------------------------------------*/
    /* Take care of indFlip itself                                            */
    /*------------------------------------------------------------------------*/
    k   = indFlip ;
    Azk = Az[k]++ ;
    if(w[k]) // wk = 1
    {
        if(Azk == 0)
        {
            IdxGroupsByC_Insert(w1Az1_groups, k) ;
        }
        else if(Azk == 1)
        {
            IdxGroupsByC_Delete(w1Az1_groups, k) ;
            IdxGroupsByC_Insert(w1Az2_groups, k) ;
        }
        else if(Azk == 2)
        {
            IdxGroupsByC_Delete(w1Az2_groups, k) ;
        }
    }
    else // wk = 0
    {
        if(Azk == 0)
        {
            IdxGroupsByC_Delete(w0Az0_groups, k) ;
            if(c[k] == gamma)
            {
                Direct_Access_List_Insert(w0Az1_list, k) ;
            }
        }
        else if(Azk == 1)
        {
            if(c[k] == gamma)
            {
                Direct_Access_List_Delete_By_Value(w0Az1_list, k) ;
            }
        }
    }

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    return ;
}


// ASSUME: that Upflip_v5 happens in z
// ASSUME: the Aw_indFlip is either 1 or 2
void Downflip_v5
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtz,
    VSPINT                       *wgtw, // used for debugging
    IdxGroupsByC         *z0Aw0_groups,
    Direct_Access_List     *z0Aw1_list,
    IdxGroupsByC         *z1Aw1_groups,
    IdxGroupsByC         *z1Aw2_groups,

    IdxGroupsByC         *w0Az0_groups,
    Direct_Access_List     *w0Az1_list,
    IdxGroupsByC         *w1Az1_groups,
    IdxGroupsByC         *w1Az2_groups,

    VSPINT                     indFlip,
    Problem               *problem
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT n = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c = problem->graph->c ;
    VSPINT *Ai   = problem->graph->Ai ;
    VSPLONG *Ap   = problem->graph->Ap ;
    VSPINT i, k, Azk ;
    VSPLONG start, end, p ;

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    z[indFlip] = 0 ;
    *wgtz -= c[indFlip] ;

    /*------------------------------------------------------------------------*/
    // Do all changes pertaining to indFlip leaving a struct
    /*------------------------------------------------------------------------*/
    VSPINT AwindFlip = Aw[indFlip] ;
    if( AwindFlip == 1 )
    {
        IdxGroupsByC_Delete(z1Aw1_groups, indFlip) ;
    }
    else if( AwindFlip == 2)
    {
        IdxGroupsByC_Delete(z1Aw2_groups, indFlip) ;
    }

    /*------------------------------------------------------------------------*/
    // Do all changes pertaining to indFlip entering a struct
    /*------------------------------------------------------------------------*/
    if(AwindFlip == 1 && c[indFlip] == gamma)
    {
        Direct_Access_List_Insert(z0Aw1_list, indFlip) ;
    }


    /*------------------------------------------------------------------------*/
    // Do all changes pertaining Awk increasing in N(indFlip) \cup {indFlip}
    /*------------------------------------------------------------------------*/
    start = Ap[indFlip] ;
    end   = Ap[indFlip + 1] ;
    for(p = start ; p < end ; p++)
    {
        k   = Ai[p] ;
        Azk = Az[k]-- ;
        if(w[k]) // wk = 1
        {
            if(Azk == 1)
            {
                IdxGroupsByC_Delete(w1Az1_groups, k) ;
            }
            else if(Azk == 2)
            {
                IdxGroupsByC_Delete(w1Az2_groups, k) ;
                IdxGroupsByC_Insert(w1Az1_groups, k) ;
            }
            else if(Azk == 3)
            {
                IdxGroupsByC_Insert(w1Az2_groups, k) ;
            }
        }
        else // wk = 0
        {
            if(Azk == 1)
            {
                IdxGroupsByC_Insert(w0Az0_groups, k) ;
                if(c[k] == gamma)
                {
                    Direct_Access_List_Delete_By_Value(w0Az1_list, k) ;
                }
            }
            else if(Azk == 2)
            {
                if(c[k] == gamma)
                {
                    Direct_Access_List_Insert(w0Az1_list, k) ;
                }
            }
        }
    }


    /*------------------------------------------------------------------------*/
    /* Take care of indFlip itself                                            */
    /*------------------------------------------------------------------------*/
    // #ifdef DEBUG
    //     IdxGroupsByC_Check(z0Aw0_groups, n) ;
    //     IdxGroupsByC_Check(z1Aw1_groups, n) ;
    //     IdxGroupsByC_Check(z1Aw2_groups, n) ;

    //     IdxGroupsByC_Check(w0Az0_groups, n) ;
    //     IdxGroupsByC_Check(w1Az1_groups, n) ;
    //     IdxGroupsByC_Check(w1Az2_groups, n) ;

    //     z0Aw1_List_Check_v5(z0Aw1_list, c, z, Aw, gamma, n) ;
    //     z0Aw1_List_Check_v5(w0Az1_list, c, w, Az, gamma, n) ;
    // #endif

    k   = indFlip ;
    Azk = Az[k]-- ;
    if(w[k]) // wk = 1
    {
        if(Azk == 1)
        {
            IdxGroupsByC_Delete(w1Az1_groups, k) ;
        }
        else if(Azk == 2)
        {
            IdxGroupsByC_Delete(w1Az2_groups, k) ;
            IdxGroupsByC_Insert(w1Az1_groups, k) ;
        }
        else if(Azk == 3)
        {
            IdxGroupsByC_Insert(w1Az2_groups, k) ;
        }
    }
    else // wk = 0
    {
        if(Azk == 1)
        {
            IdxGroupsByC_Insert(w0Az0_groups, k) ;
            if(c[k] == gamma)
            {
                Direct_Access_List_Delete_By_Value(w0Az1_list, k) ;
            }
        }
        else if(Azk == 2)
        {
            if(c[k] == gamma)
            {
                Direct_Access_List_Insert(w0Az1_list, k) ;
            }
        }
    }

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif
    return ;
}


void Initialize_All_Structs_v5
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    IdxGroupsByC         *x0Ay0_groups,
    Direct_Access_List     *x0Ay1_list,
    IdxGroupsByC         *x1Ay1_groups,
    IdxGroupsByC         *x1Ay2_groups,

    IdxGroupsByC         *y0Ax0_groups,
    Direct_Access_List     *y0Ax1_list,
    IdxGroupsByC         *y1Ax1_groups,
    IdxGroupsByC         *y1Ax2_groups,

    Problem               *problem
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT i, xi, yi, Axi, Ayi, ci ;
    VSPINT *c    = problem->graph->c ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    for(i = 0 ; i < n ; i++)
    {
        xi = x[i] ;
        yi = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;
        ci = c[i] ;

        if(xi)
        {
            if(Ayi == 1)
            {
                IdxGroupsByC_Insert(x1Ay1_groups, i) ;
            }
            else if(Ayi == 2)
            {
                IdxGroupsByC_Insert(x1Ay2_groups, i) ;
            }
        }
        else // xi = 0
        {
            if(Ayi == 0)
            {
                IdxGroupsByC_Insert(x0Ay0_groups, i) ;
            }
            else if(Ayi == 1 && ci == gamma)
            {
                Direct_Access_List_Insert(x0Ay1_list, i) ;
            }
        }

        if(yi)
        {
            if(Axi == 1)
            {
                IdxGroupsByC_Insert(y1Ax1_groups, i) ;
            }
            else if(Axi == 2)
            {
                IdxGroupsByC_Insert(y1Ax2_groups, i) ;
            }
        }
        else // yi = 0
        {
            if(Axi == 0)
            {
                IdxGroupsByC_Insert(y0Ax0_groups, i) ;
            }
            else if(Axi == 1 && ci == gamma)
            {
                Direct_Access_List_Insert(y0Ax1_list, i) ;
            }
        }
    }

    #ifdef DEBUG
        Check_All_Data_v5
        (
            x, Ax, wgtx, y, Ay, wgty,
            x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
            y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
            problem
        ) ;
    #endif

    return ;
}



void Check_All_Data_v5
(
    VSPINT                          *x,
    VSPINT                         *Ax,
    VSPINT                        wgtx,
    VSPINT                          *y,
    VSPINT                         *Ay,
    VSPINT                        wgty,
    IdxGroupsByC         *x0Ay0_groups,
    Direct_Access_List     *x0Ay1_list,
    IdxGroupsByC         *x1Ay1_groups,
    IdxGroupsByC         *x1Ay2_groups,

    IdxGroupsByC         *y0Ax0_groups,
    Direct_Access_List     *y0Ax1_list,
    IdxGroupsByC         *y1Ax1_groups,
    IdxGroupsByC         *y1Ax2_groups,

    Problem               *problem
)
{
    SET_FUNCTION_COUNTER ;
    WARN_IF_NOT_DEBUG_MODE ;

    #ifdef DEBUG
    VSPINT i ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;

    /*------------------------------------------------------------------------*/
    // Load the data into states and use check_states
    /*------------------------------------------------------------------------*/
    Point xPoint, yPoint ;

    xPoint.int_part = x ;
    xPoint.ind_frac = EMPTY ;
    xPoint.val_frac = 0 ;
    xPoint.n = n ;
    
    yPoint.int_part = y ;
    yPoint.ind_frac = EMPTY ;
    yPoint.val_frac = 0 ;
    yPoint.n = n ;
    
    Az Axstruct, Aystruct ;

    Axstruct.int_part    = Ax ;
    Axstruct.ind_frac_z  = EMPTY ;
    Axstruct.val_frac_z  = 0 ;
    Axstruct.n           = n ;
    Axstruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Axstruct.one         = 1 ;


    Aystruct.int_part    = Ay ;
    Aystruct.ind_frac_z  = EMPTY ;
    Aystruct.val_frac_z  = 0 ;
    Aystruct.n           = n ;
    Aystruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Aystruct.one         = 1 ;

    State X, Y ;

    X.point  = &xPoint ;
    X.Az     = &Axstruct ;
    X.weight = wgtx ;

    Y.point  = &yPoint ;
    Y.Az     = &Aystruct ;
    Y.weight = wgty ;

    check_states(&X, &Y, problem) ;


    if(wgtx < 0 || wgty < 0)
    {
        ERR("Negative weights") ;
    }


    /*------------------------------------------------------------------------*/
    // Check all the groups and direct access lists
    /*------------------------------------------------------------------------*/
    // Check x groups
    IdxGroupsByC_Check(x0Ay0_groups, n) ;
    IdxGroupsByC_Check(x1Ay1_groups, n) ;
    IdxGroupsByC_Check(x1Ay2_groups, n) ;

    // Check y groups
    IdxGroupsByC_Check(y0Ax0_groups, n) ;
    IdxGroupsByC_Check(y1Ax1_groups, n) ;
    IdxGroupsByC_Check(y1Ax2_groups, n) ;

    // Check x lists
    Direct_Access_List_Check(x0Ay1_list) ;
    z0Aw1_List_Check_v5(x0Ay1_list, c, x, Ay, gamma, n) ;

    // Check y lists
    Direct_Access_List_Check(y0Ax1_list) ;
    z0Aw1_List_Check_v5(y0Ax1_list, c, y, Ax, gamma, n) ;

    /*------------------------------------------------------------------------*/
    // Check that there are no i with xi = 1 and Ayi >= 3 or 
    // yi = 1 and Axi >= 3
    /*------------------------------------------------------------------------*/
    for(i = 0 ; i < n ; i++)
    {
        if(x[i] == 1 && Ay[i] >= 3)
        {
            ERR("Downflip_v5 in x with penalty >= 3") ;
        }

        if(y[i] == 1 && Ax[i] >= 3)
        {
            ERR("Downflip_v5 in y with penalty >= 3") ;
        }
    }

    // Free memory
    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;
    #endif
    
    return ;
}

// Check that the given flip is indeed the best
// void Check_Given_Flip_Is_Optimal_v5
// (
//     VSPINT                  *x,
//     VSPINT                 *Ax,
//     VSPINT                wgtx,
//     VSPINT                  *y,
//     VSPINT                 *Ay,
//     VSPINT                wgty,
//     char          givenVarFlip,
//     VSPINT        givenIndFlip,
//     Problem       *problem
// )
// {
//     static long long unsigned count = 0 ;
//     count++ ;
//     WARN_IF_NOT_DEBUG_MODE ;

//     VSPINT i, xi, yi, Axi, Ayi, ci ;
//     VSPINT n = problem->graph->n ;
//     VSPINT gamma = problem->graph->max_vtx_cost ;
//     VSPINT *c = problem->graph->c ;
//     VSPINT  u = problem->u ;

//     Flip trueFlip ;
//     VSPINT gain, maxFeasGain = 0 ;
//     for(i = 0 ; i < n ; i++)
//     {
//         xi  = x[i] ;
//         yi  = y[i] ;
//         Axi = Ax[i] ;
//         Ayi = Ay[i] ;
//         ci  = c[i] ;
//         if(xi)
//         {
//             gain = GAIN_FLIP(DOWN, ci, gamma, Ayi) ;
//             if(gain > maxFeasGain)
//         }
//         else
//         {

//         }
//     }












//     ////////////////////////////////////////////////////////////////


//     /*------------------------------------------------------------------------*/
//     // Get the p values
//     /*------------------------------------------------------------------------*/
//     VSPINT px0 = n + 1, py0 = n + 1, px1 = -1, py1 = -1 ;

//     for(i = 0 ; i < n ; i++)
//     {
//         xi = x[i] ;
//         yi = y[i] ;
//         Axi = Ax[i] ;
//         Ayi = Ay[i] ;

//         // Record px's
//         if(xi)
//         {
//             if(Ayi > px1)
//             {
//                 px1 = Ayi ;
//             }
//         }
//         else
//         {
//             if(Ayi < px0)
//             {
//                 px0 = Ayi ;
//             }
//         }

//         // Record py's
//         if(yi)
//         {
//             if(Axi > py1)
//             {
//                 py1 = Axi ;
//             }
//         }
//         else
//         {
//             if(Axi < py0)
//             {
//                 py0 = Axi ;
//             }
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Record the c's and an index for each c
//     // Also record:
//     // x0witness, y0witness, x1witness, y1witness - indices achieving the
//     // optimal p and c values
//     /*------------------------------------------------------------------------*/
//     VSPINT cx0 = -1, cy0 = -1, cx1 = gamma + 1, cy1 = gamma + 1 ;
//     VSPINT x0witness = EMPTY, y0witness = EMPTY ;
//     VSPINT x1witness = EMPTY, y1witness = EMPTY ;

//     for(i = 0 ; i < n ; i++)
//     {
//         xi = x[i] ;
//         yi = y[i] ;
//         Axi = Ax[i] ;
//         Ayi = Ay[i] ;
//         ci = c[i] ;

//         // Deal with x
//         if(xi == 0 && Ayi == px0)
//         {
//             if(ci > cx0 )
//             {
//                 cx0 = ci ;
//                 x0witness = i ;
//             }
//         }

//         if(xi == 1 && Ayi == px1)
//         {
//             if(ci < cx1 )
//             {
//                 cx1 = ci ;
//                 x1witness = i ;
//             }
//         }

//         // Deal with y
//         if(yi == 0 && Axi == py0)
//         {
//             if(ci > cy0 )
//             {
//                 cy0 = ci ;
//                 y0witness = i ;
//             }
//         }

//         if(yi == 1 && Axi == py1)
//         {
//             if(ci < cy1 )
//             {
//                 cy1 = ci ;
//                 y1witness = i ;
//             }
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Go through the 4 possible flips and record the maximum feasible gain
//     /*------------------------------------------------------------------------*/
//     VSPINT trueBestGain    = 0 ;
//     VSPINT trueBestGainInd = EMPTY ;
//     char   trueBestGainVar = '\0' ;
//     VSPINT trueBestGainDir = 0 ; // 1 for Upflip_v5, -1 for Downflip_v5

//     VSPINT x0gain = cx0 - gamma * px0 ;
//     VSPINT y0gain = cy0 - gamma * py0 ;
//     VSPINT x1gain = gamma * px1 - cx1 ;
//     VSPINT y1gain = gamma * py1 - cy1 ;

//     if( FEASIBLE_UPFLIP(wgtx, cx0, u) )
//     {
//         if(x0gain > trueBestGain)
//         {
//             trueBestGain    = x0gain ;
//             trueBestGainInd = x0witness ;
//             trueBestGainVar = 'x' ;
//             trueBestGainDir = 1 ;
//         }
//     }

//     if( FEASIBLE_UPFLIP(wgty, cy0, u) )
//     {
//         if(y0gain > trueBestGain)
//         {
//             trueBestGain    = y0gain ;
//             trueBestGainInd = y0witness ;
//             trueBestGainVar = 'y' ;
//             trueBestGainDir = 1 ;
//         }
//     }

//     if( FEASIBLE_DNFLIP(wgtx, cx1) )
//     {
//         if(x1gain > trueBestGain)
//         {
//             trueBestGain    = x1gain ;
//             trueBestGainInd = x1witness ;
//             trueBestGainVar = 'x' ;
//             trueBestGainDir = -1 ;
//         }
//     }

//     if( FEASIBLE_DNFLIP(wgty, cy1) )
//     {
//         if(y1gain > trueBestGain)
//         {
//             trueBestGain    = y1gain ;
//             trueBestGainInd = y1witness ;
//             trueBestGainVar = 'y' ;
//             trueBestGainDir = -1 ;
//         }
//     }

//     /*------------------------------------------------------------------------*/
//     // Compute the gain of the given flip. Load the given and best flips into
//     // flip structs for easier debugging
//     /*------------------------------------------------------------------------*/

//     // Compute the gain of the given var
//     VSPINT givenGain = EMPTY ;
//     VSPINT dir ;

//     if(givenIndFlip != EMPTY)
//     {
//         if(givenVarFlip == 'x')
//         {
//             dir = (1 - 2 * x[givenIndFlip]) ;
//             givenGain = dir * (c[givenIndFlip] - gamma *(Ay[givenIndFlip] )) ;
//         }
//         else if(givenVarFlip == 'y')
//         {
//             dir = (1 - 2 * y[givenIndFlip]) ;
//             givenGain = dir * (c[givenIndFlip] - gamma *(Ax[givenIndFlip] )) ;
//         }
//     }

//     // Load the given flip into a Flip struct for easier debugging
//     Flip givenFlip = {0} ;
//     givenFlip.gamma = gamma ;
//     givenFlip.i    = givenIndFlip ;
//     givenFlip.var  = givenVarFlip ;
//     givenFlip.gain = givenGain ;
//     givenFlip.u = u ;
//     if(givenVarFlip == 'x')
//     {
//         VSP_ASSERT(givenIndFlip != EMPTY) ;
//         givenFlip.wgt = wgtx ;
//         givenFlip.zi   = x[givenIndFlip] ;
//         givenFlip.Awi  = Ay[givenIndFlip] ;
//         givenFlip.ci   = c[givenIndFlip] ;
//         if(x[givenIndFlip] == 0)
//         {
//             givenFlip.dir  = 1 ;
//         }
//         else
//         {
//             givenFlip.dir  = -1 ;
//         }
//     }
//     else if(givenVarFlip == 'y')
//     {
//         VSP_ASSERT(givenIndFlip != EMPTY) ;
//         givenFlip.wgt = wgty ;
//         givenFlip.zi   = y[givenIndFlip] ;
//         givenFlip.Awi  = Ax[givenIndFlip] ;
//         givenFlip.dir  = 1 - 2 * y[givenIndFlip] ;
//         givenFlip.ci   = c[givenIndFlip] ;
//         if(y[givenIndFlip] == 0)
//         {
//             givenFlip.dir  = 1 ;
//         }
//         else
//         {
//             givenFlip.dir  = -1 ;
//         }
//     }

//     // Load the true optimal flip into its own flip struct for easier debugging
//     Flip trueFlip = {0} ;

//     trueFlip.gamma = gamma ;
//     trueFlip.i = trueBestGainInd ;
//     trueFlip.gain = trueBestGain ;
//     trueFlip.var = trueBestGainVar ;
//     if(trueBestGainVar == 'x') trueFlip.wgt = wgtx ;
//     else if(trueBestGainVar == 'y') trueFlip.wgt = wgty ;
//     trueFlip.u = u ;

//     if(trueBestGainInd != EMPTY)
//     {
//         trueFlip.ci = c[trueBestGainInd] ;
        
//         if(trueBestGainVar == 'x')
//         {
//             trueFlip.zi  = x[trueBestGainInd] ;
//             trueFlip.Awi = Ay[trueBestGainInd] ;
//             trueFlip.dir = ( (1 - 2 * x[trueBestGainInd]) == 1 ) ? 1 : -1 ;
//         }
//         else
//         {
//             trueFlip.zi  = y[trueBestGainInd] ;
//             trueFlip.Awi = Ax[trueBestGainInd] ;
//             trueFlip.dir = ( (1 - 2 * y[trueBestGainInd]) == 1 ) ? 1 : -1 ;
//         }
//     }


//     /*------------------------------------------------------------------------*/
//     // Check the given flip against the true best flip
//     /*------------------------------------------------------------------------*/
//     // For debugging:
//     // trueFlip struct
//     // givenFlip struct

//     VSPINT givenFlipExists  = (givenIndFlip != EMPTY) ;
//     VSPINT trueBestExists   = (trueBestGainInd != EMPTY) ;

//     // Either both must be true or both must be false simultaneously
//     if( (givenFlipExists + trueBestExists != 0) && (givenFlipExists + trueBestExists != 2) ) 
//     {
//         ERR("Mismatch between given flip and true best flip") ;
//     }

//     if(trueBestExists)
//     {
//         if(givenFlip.gain != trueFlip.gain)
//         {
//             ERR("Mismatch between gain of given flip and true best flip") ;
//         }
//     }

//     return ;
// }


void Check_Inputs_v5
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    Problem       *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    
    WARN_IF_NOT_DEBUG_MODE ;

    VSPINT i ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;

    /*------------------------------------------------------------------------*/
    // Load the data into states and use check_states
    /*------------------------------------------------------------------------*/
    Point xPoint, yPoint ;

    xPoint.int_part = x ;
    xPoint.ind_frac = EMPTY ;
    xPoint.val_frac = 0 ;
    xPoint.n = n ;
    
    yPoint.int_part = y ;
    yPoint.ind_frac = EMPTY ;
    yPoint.val_frac = 0 ;
    yPoint.n = n ;
    
    Az Axstruct, Aystruct ;

    Axstruct.int_part    = Ax ;
    Axstruct.ind_frac_z  = EMPTY ;
    Axstruct.val_frac_z  = 0 ;
    Axstruct.n           = n ;
    Axstruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Axstruct.one         = 1 ;


    Aystruct.int_part    = Ay ;
    Aystruct.ind_frac_z  = EMPTY ;
    Aystruct.val_frac_z  = 0 ;
    Aystruct.n           = n ;
    Aystruct.col_frac_z  = calloc(n, sizeof(VSPINT)) ;
    Aystruct.one         = 1 ;

    State X, Y ;

    X.point  = &xPoint ;
    X.Az     = &Axstruct ;
    X.weight = wgtx ;

    Y.point  = &yPoint ;
    Y.Az     = &Aystruct ;
    Y.weight = wgty ;

    check_states(&X, &Y, problem) ;

    if(wgtx < 0 || wgty < 0)
    {
        ERR("Negative weights") ;
    }

    /*------------------------------------------------------------------------*/
    // Check that there are no i with xi = 1 and Ayi >= 3 or 
    // yi = 1 and Axi >= 3
    /*------------------------------------------------------------------------*/
    for(i = 0 ; i < n ; i++)
    {
        if(x[i] == 1 && Ay[i] >= 3)
        {
            ERR("Downflip_v5 in x with penalty >= 3") ;
        }

        if(y[i] == 1 && Ax[i] >= 3)
        {
            ERR("Downflip_v5 in y with penalty >= 3") ;
        }
    }

    // Free memory
    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;

    return ;
}


// In addition to checking the data,check that there are no strictly beneficial
// flips 
void Check_Outputs_v5
(
    VSPINT                  *x,
    VSPINT                 *Ax,
    VSPINT                wgtx,
    VSPINT                  *y,
    VSPINT                 *Ay,
    VSPINT                wgty,
    Problem       *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;
    Check_Inputs_v5(x, Ax, wgtx, y, Ay, wgty, problem) ;

    VSPINT n = problem->graph->n ;
    VSPINT u = problem->u ;
    VSPINT *c = problem->graph->c ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    VSPINT ps[4], cs[4], witness[4], isFeas[4], gains[4] ;
    Compute_ps_And_cs_Dbg_v5(x, Ax, wgtx, y, Ay, wgty, problem, ps, cs, witness) ;

    Compute_Feasibility_And_Gains_v5
    (
        Ax, Ay, wgtx, wgty, problem,
        witness, isFeas, gains
    ) ;

    for(VSPINT k = 0 ; k < 4 ; k++)
    {
        if(isFeas[k] && gains[k] > 0)
        {
            ERR("Quit strictly beneficial flips with one still available") ;
        }
    }

    // WARNING: this may give false positive, as we exit when structs tell
    // us there are no more flips, but remember non-root heap node may
    // correspond to a beneficial Upflip_v5 that the structs can't see
    // for(VSPINT i = 0 ; i < n ; i++)
    // {
    //     if(x[i])
    //     {
    //         if(GAIN_FLIP(-1, c[i], gamma, Ay[i]) > 0 && FEASIBLE_DNFLIP(wgtx, c[i]))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Downflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    //     else
    //     {
    //         if(GAIN_FLIP(1, c[i], gamma, Ay[i]) > 0 && FEASIBLE_UPFLIP(wgtx, c[i], u))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Upflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }

    //     if(y[i])
    //     {
    //         if(GAIN_FLIP(-1, c[i], gamma, Ax[i]) > 0 && FEASIBLE_DNFLIP(wgty, c[i]))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Downflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    //     else
    //     {
    //         if(GAIN_FLIP(1, c[i], gamma, Ax[i]) > 0 && FEASIBLE_UPFLIP(wgty, c[i], u))
    //         {
    //             // If I get this error, it may be a false positive. Read the 
    //             // warning above
    //             ERR("Feasible strictly beneficial Upflip_v5 in y at end of flips") ;
    //             break ;
    //         }
    //     }
    // }
}



// ASSUMED: No assumptions on maxGroupSizes, it is set to 0 at beginning. It
// may start out as garbage
void Compute_Max_Group_Sizes_v5
(
    VSPINT *maxGroupSizes,
    VSPINT *c,
    VSPINT gamma,
    VSPINT n
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i ;
    for(i = 0 ; i <= gamma ; i++)
    {
        maxGroupSizes[i] = 0 ;
    }

    for(i = 0 ; i < n ; i++)
    {
        VSP_ASSERT( (IN_RANGE(c[i], 1, gamma + 1)) ) ;
        maxGroupSizes[ c[i] ]++ ;
    }
    return ;
}

void Do_Flip_v5
(
    VSPINT                          *z,
    VSPINT                         *Az,
    VSPINT                       *wgtz,
    VSPINT                          *w,
    VSPINT                         *Aw,
    VSPINT                       *wgtw, // used for debugging only
    IdxGroupsByC         *z0Aw0_groups,
    Direct_Access_List     *z0Aw1_list,
    IdxGroupsByC         *z1Aw1_groups,
    IdxGroupsByC         *z1Aw2_groups,

    IdxGroupsByC         *w0Az0_groups,
    Direct_Access_List     *w0Az1_list,
    IdxGroupsByC         *w1Az1_groups,
    IdxGroupsByC         *w1Az2_groups,

    VSPINT                     indFlip,
    Problem               *problem
)
{
    SET_FUNCTION_COUNTER ;

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    if(z[indFlip])
    {
        Downflip_v5
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

            indFlip, problem
        ) ;
    }
    else
    {
        Upflip_v5
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

            indFlip, problem
        ) ;
    }

    #ifdef DEBUG
        Check_All_Data_v5
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    return ;
}



// Assign the lighter variable to z and the heavier one to w
void Assign_By_Weight_v5
(
    VSPINT **z,
    VSPINT **w,
    VSPINT **Az,
    VSPINT **Aw,
    VSPINT **wgtz,
    VSPINT **wgtw,

    IdxGroupsByC **z0Aw0_groups,
    IdxGroupsByC **z1Aw1_groups,
    IdxGroupsByC **z1Aw2_groups,
    Direct_Access_List **z0Aw1_list,

    IdxGroupsByC **w0Az0_groups,
    IdxGroupsByC **w1Az1_groups,
    IdxGroupsByC **w1Az2_groups,
    Direct_Access_List **w0Az1_list,

    VSPINT *x,
    VSPINT *y,
    VSPINT *Ax,
    VSPINT *Ay,
    VSPINT *wgtx,
    VSPINT *wgty,

    IdxGroupsByC *x0Ay0_groups,
    IdxGroupsByC *x1Ay1_groups,
    IdxGroupsByC *x1Ay2_groups,
    Direct_Access_List *x0Ay1_list,

    IdxGroupsByC *y0Ax0_groups,
    IdxGroupsByC *y1Ax1_groups,
    IdxGroupsByC *y1Ax2_groups,
    Direct_Access_List *y0Ax1_list
)
{
    if (*wgtx < *wgty) {
        *z = x;
        *w = y;
        *Az = Ax;
        *Aw = Ay;
        *wgtz = wgtx;
        *wgtw = wgty;

        *z0Aw0_groups = x0Ay0_groups;
        *z1Aw1_groups = x1Ay1_groups;
        *z1Aw2_groups = x1Ay2_groups;
        *z0Aw1_list = x0Ay1_list;

        *w0Az0_groups = y0Ax0_groups;
        *w1Az1_groups = y1Ax1_groups;
        *w1Az2_groups = y1Ax2_groups;
        *w0Az1_list = y0Ax1_list;
    }
    else
    {
        *z = y;
        *w = x;
        *Az = Ay;
        *Aw = Ax;
        *wgtz = wgty;
        *wgtw = wgtx;

        *z0Aw0_groups = y0Ax0_groups;
        *z1Aw1_groups = y1Ax1_groups;
        *z1Aw2_groups = y1Ax2_groups;
        *z0Aw1_list = y0Ax1_list;

        *w0Az0_groups = x0Ay0_groups;
        *w1Az1_groups = x1Ay1_groups;
        *w1Az2_groups = x1Ay2_groups;
        *w0Az1_list = x0Ay1_list;
    }
}
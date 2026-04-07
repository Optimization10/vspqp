#include "core.h"
#include "greedy_flips_internal.h"
#include "utils.h"

void z0Aw1_List_Check
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

// Assume that upflip happens in z
void Upflip
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
        Check_All_Data
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
        Check_All_Data
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif

    return ;
}


// ASSUME: that upflip happens in z
// ASSUME: the Aw_indFlip is either 1 or 2
void Downflip
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
        Check_All_Data
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
        Check_All_Data
        (
            z, Az, *wgtz, w, Aw, *wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,
            problem
        ) ;
    #endif
    return ;
}


void Initialize_All_Structs
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
        Check_All_Data
        (
            x, Ax, wgtx, y, Ay, wgty,
            x0Ay0_groups, x0Ay1_list, x1Ay1_groups, x1Ay2_groups,
            y0Ax0_groups, y0Ax1_list, y1Ax1_groups, y1Ax2_groups,
            problem
        ) ;
    #endif

    return ;
}



void Check_All_Data
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
    z0Aw1_List_Check(x0Ay1_list, c, x, Ay, gamma, n) ;

    // Check y lists
    Direct_Access_List_Check(y0Ax1_list) ;
    z0Aw1_List_Check(y0Ax1_list, c, y, Ax, gamma, n) ;

    /*------------------------------------------------------------------------*/
    // Check that there are no i with xi = 1 and Ayi >= 3 or 
    // yi = 1 and Axi >= 3
    /*------------------------------------------------------------------------*/
    for(i = 0 ; i < n ; i++)
    {
        if(x[i] == 1 && Ay[i] >= 3)
        {
            ERR("Downflip in x with penalty >= 3") ;
        }

        if(y[i] == 1 && Ax[i] >= 3)
        {
            ERR("Downflip in y with penalty >= 3") ;
        }
    }

    // Free memory
    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;
    #endif
    
    return ;
}

void Check_Inputs
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
            ERR("Downflip in x with penalty >= 3") ;
        }

        if(y[i] == 1 && Ax[i] >= 3)
        {
            ERR("Downflip in y with penalty >= 3") ;
        }
    }

    // Free memory
    free(Axstruct.col_frac_z) ;
    free(Aystruct.col_frac_z) ;

    return ;
}


// In addition to checking the data,check that there are no strictly beneficial
// flips 
void Check_Outputs
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
    Check_Inputs(x, Ax, wgtx, y, Ay, wgty, problem) ;

    VSPINT n = problem->graph->n ;
    VSPINT u = problem->u ;
    VSPINT *c = problem->graph->c ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    VSPINT ps[4], cs[4], witness[4], isFeas[4], gains[4] ;
    Compute_ps_And_cs_Dbg(x, Ax, wgtx, NULL, y, Ay, wgty, NULL, problem, ps, cs, witness) ;

    Compute_Feasibility_And_Gains
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
    Problem           *problem
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;
    Check_Inputs(x, Ax, wgtx, y, Ay, wgty, problem) ;

    VSPINT n = problem->graph->n ;
    VSPINT u = problem->u ;
    VSPINT *c = problem->graph->c ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    VSPINT ps[4], cs[4], witness[4], isFeas[4], gains[4] ;
    Compute_ps_And_cs_Dbg(x, Ax, wgtx, NULL, y, Ay, wgty, NULL, problem, ps, cs, witness) ;

    Compute_Feasibility_And_Gains
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
}

// ASSUMED: No assumptions on maxGroupSizes, it is set to 0 at beginning. It
// may start out as garbage
void Compute_Max_Group_Sizes
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

void Do_Flip
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
    Problem                   *problem
)
{
    SET_FUNCTION_COUNTER ;
    if(z[indFlip])
    {
        Downflip
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

            indFlip, problem
        ) ;
    }
    else
    {
        Upflip
        (
            z, Az, w, Aw, wgtz, wgtw,
            z0Aw0_groups, z0Aw1_list, z1Aw1_groups, z1Aw2_groups,
            w0Az0_groups, w0Az1_list, w1Az1_groups, w1Az2_groups,

            indFlip, problem
        ) ;
    }
    return ;
}






void Compute_ps_And_cs_Dbg
(
    VSPINT                    *x,
    VSPINT                   *Ax,
    VSPINT                  wgtx,
    VSPINT              *xlocked, // Ignored if NULL
    VSPINT                    *y,
    VSPINT                   *Ay,
    VSPINT                  wgty,
    VSPINT              *ylocked, // Ignored if NULL
    Problem             *problem,
    VSPINT                 ps[4],
    VSPINT                 cs[4],
    VSPINT          witnesses[4]
)
{
    static long long unsigned count = 0 ;
    count++ ;

    WARN_IF_NOT_DEBUG_MODE ;

    // Check that either one or both variables are locked
    if((xlocked == NULL && ylocked != NULL) ||
       (xlocked != NULL && ylocked == NULL)    )
    {
        ERR("Only one variable locked") ;
    }

    VSPINT i, xi, yi, Axi, Ayi, ci ;
    VSPINT n     = problem->graph->n ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPINT *c    = problem->graph->c ;
    VSPINT  u    = problem->u ;

    /*------------------------------------------------------------------------*/
    // Get the p values
    /*------------------------------------------------------------------------*/
    VSPINT px0 = n + 1, py0 = n + 1, px1 = -1, py1 = -1 ;

    for(i = 0 ; i < n ; i++)
    {
        xi  = x[i] ;
        yi  = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;

        if( UNLOCKED(xlocked, i) )
        {
            // Record px's
            if(xi)
            {
                if(Ayi > px1)
                {
                    px1 = Ayi ;
                }
            }
            else
            {
                if(Ayi < px0)
                {
                    px0 = Ayi ;
                }
            }
        }

        if( UNLOCKED(ylocked, i) )
        {
            // Record py's
            if(yi)
            {
                if(Axi > py1)
                {
                    py1 = Axi ;
                }
            }
            else
            {
                if(Axi < py0)
                {
                    py0 = Axi ;
                }
            }
        }
    }

    /*------------------------------------------------------------------------*/
    // Record the c's and an index for each c
    // Also record:
    // x0witness, y0witness, x1witness, y1witness - indices achieving the
    // optimal p and c values
    /*------------------------------------------------------------------------*/
    VSPINT cx0 = 0, cy0 = 0, cx1 = gamma + 1, cy1 = gamma + 1 ;
    VSPINT x0witness = EMPTY, y0witness = EMPTY ;
    VSPINT x1witness = EMPTY, y1witness = EMPTY ;

    for(i = 0 ; i < n ; i++)
    {
        xi = x[i] ;
        yi = y[i] ;
        Axi = Ax[i] ;
        Ayi = Ay[i] ;
        ci = c[i] ;

        // Deal with x
        if(xi == 0 && Ayi == px0 && UNLOCKED(xlocked, i) )
        {
            if(ci > cx0 )
            {
                cx0 = ci ;
                x0witness = i ;
            }
        }

        if(xi == 1 && Ayi == px1 && UNLOCKED(xlocked, i))
        {
            if(ci < cx1 )
            {
                cx1 = ci ;
                x1witness = i ;
            }
        }

        // Deal with y
        if(yi == 0 && Axi == py0 && UNLOCKED(ylocked, i))
        {
            if(ci > cy0 )
            {
                cy0 = ci ;
                y0witness = i ;
            }
        }

        if(yi == 1 && Axi == py1 && UNLOCKED(ylocked, i) )
        {
            if(ci < cy1 )
            {
                cy1 = ci ;
                y1witness = i ;
            }
        }
    }

    // Prepare return
    ps[X0] = px0 ;
    ps[X1] = px1 ;
    ps[Y0] = py0 ;
    ps[Y1] = py1 ;

    cs[X0] = cx0 ;
    cs[X1] = cx1 ;
    cs[Y0] = cy0 ;
    cs[Y1] = cy1 ;

    witnesses[X0] = x0witness ;
    witnesses[X1] = x1witness ;
    witnesses[Y0] = y0witness ;
    witnesses[Y1] = y1witness ;

    return ;
}




/*------------------------------------------------------------------------------
Record the gain even if it's negative

inds[k] - k is for the encoding, inds[k] is a vertex number. inds[k] may be
          EMPTY

retIsFeas[k] - if inds[k] != EMPTY and the flip is feasible, the it's 1
               otherwise, it's 0

retGains[k] - if retIsFeas[k] = 1, then this is the gain of the flip. If
              retIsFeas[k] = 0, then this is garbage
------------------------------------------------------------------------------*/
void Compute_Feasibility_And_Gains
(
    VSPINT                *Ax,
    VSPINT                *Ay,
    VSPINT               wgtx,
    VSPINT               wgty,
    Problem      *problem,
    VSPINT             inds[4],
    VSPINT        retIsFeas[4],
    VSPINT         retGains[4]
)
{
    SET_FUNCTION_COUNTER ;
    VSPINT i ;
    VSPINT n  = problem->graph->n ;
    VSPINT *c = problem->graph->c ;
    VSPINT u  = problem->u ;
    VSPINT gamma = problem->graph->max_vtx_cost ;

    for(i = 0 ; i < 4 ; i++) retIsFeas[i] = 0 ;

    if(inds[X0] != EMPTY && FEASIBLE_UPFLIP(wgtx, c[inds[X0]], u) )
    {
        retIsFeas[X0] = 1 ;
        retGains [X0] = GAIN_FLIP(1,  c[inds[X0]], gamma, Ay[inds[X0]]) ;
    }

    if(inds[Y0] != EMPTY && FEASIBLE_UPFLIP(wgty, c[inds[Y0]], u) )
    {
        retIsFeas[Y0] = 1 ;
        retGains [Y0] = GAIN_FLIP(1,  c[inds[Y0]], gamma, Ax[inds[Y0]]) ;
    }

    if(inds[X1] != EMPTY)
    {
        retIsFeas[X1] = 1 ;
        retGains [X1] = GAIN_FLIP(-1,  c[inds[X1]], gamma, Ay[inds[X1]]) ;
    }

    if(inds[Y1] != EMPTY)
    {
        retIsFeas[Y1] = 1 ;
        retGains [Y1] = GAIN_FLIP(-1,  c[inds[Y1]], gamma, Ax[inds[Y1]]) ;
    }

    return ;
}


void Do_Flip_Basic
(
    VSPINT    *z,
    VSPINT   *Az,
    VSPINT     i,
    VSPINT *wgtz,
    VSPINT    *c,
    VSPINT   *Ai,
    VSPLONG *Ap
)
{
    VSPINT dir, j ;
    VSPLONG end, p ;

    // Do the flip and update weight
    dir    = 1 - 2 * z[i] ;
    z[i]  += dir ;
    *wgtz += dir * c[i] ;

    // Update Az
    end   = Ap[i+1] ;
    for(p = Ap[i] ; p < end ; p++)
    {
        j = Ai[p] ;
        Az[ j ] += dir ;
    }
    Az[i] += dir ;
}
#include "shared_api.h"



void improve
(
    State         *X,
    State         *Y,
    Problem *problem
) ;

/* ============================================================================
   ================ relaxed_to_vertex_separator ===========================
   ============================================================================

 PURPOSE: Given binary points x, y which satisfy the balance constraints,
 convert them to binary points xbar, ybar which satisfy the balance constraints,
 the constraint xbar * A * ybar = 0, and have an objective value at least 
 as good as x, y

 INPUTS:
 a) A binary point x satisfying the balance constraint
 b) A binary point y satisfying the balance constraint

 OUTPUTS:
 a) Binary points xbar, ybar which satisfy the balance constraints, the
    constraint xbar * A * ybar = 0, and have an objective value at least 
    as good as x, y
*/

void relaxed_to_vertex_separator
(
    State            *X,
    State            *Y,
    Problem    *problem
)
{
    // printf("\nStarting Relaxed To VSP") ;
    fflush(stdout) ;
    #ifdef DEBUG
        check_states(X, Y, problem) ;
    #endif


    TIMER_START ;

    /* Local vars */
    VSPINT *x, *y, *Ax, *Ay, *c, *ind_bad_x, *ind_bad_y ;
    VSPINT i, j, n, xi, yi, xAy, sum_x, sum_y, num_bad_x, num_bad_y;
    xAy = 0 ;
    sum_x = 0 ;
    sum_y = 0 ;
    n = problem->Graph->n ;
    Ax = X->Az->int_part ;
    Ay = Y->Az->int_part ;
    x = X->point->int_part ;
    y = Y->point->int_part ;
    c = problem->Graph->c ;

    /* Set the fractional part to 0 or 1 */
    Force_To_Binary_Point(X, Y, problem) ;

    #ifdef DEBUG
        static unsigned long long count = 0 ;
        count += 1 ;

        /* Check the inputs */
        check_states(X, Y, problem) ;

        /* Compute the initial objective */
        VSPFLOAT objin = compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        VSPINT obj_old = 0 ;
        VSPINT gamma = problem->Graph->max_vtx_cost ;

        for (i = 0 ; i < n ; i++) {
            obj_old += c [i] * ( x [i] + y [i] ) - gamma * x [i] * Ay [i] ;
        }
    #endif

    /* ind_bad_x is array of indices i of x with x_i = 1 and (Ay)_i >= 1 */
    ind_bad_x = malloc(n * sizeof(VSPINT)) ;
    ind_bad_y = malloc(n * sizeof(VSPINT)) ;
    num_bad_x = 0 ;
    num_bad_y = 0 ;

    /* Make a preliminary pass to record all the info */
    for (i = 0 ; i < n ; i++) {
        xi = x [i] ;
        yi = y [i] ;

        xAy += xi * Ay [i] ;

        if (xi) {
            sum_x += 1 ;
            if ( Ay [i] ) {
                ind_bad_x [num_bad_x] = i;
                num_bad_x += 1 ;
            }
        }
        if (yi) {
            sum_y += 1 ;
            if ( Ax [i] ) {
                ind_bad_y [num_bad_y] = i;
                num_bad_y += 1 ;
            }
        }
    }

    bool terribleSituation = false ;
    while (xAy > 0) {
        if (sum_x > 1) {
            j = ind_bad_x [num_bad_x - 1] ;

            /* Set xj to 0. Update Ax and xAy accordingly */
            x [j] = 0 ;
            num_bad_x -= 1 ;
            sum_x-- ;

            Az_pm_col_of_A(Ax, j, MINUS, problem) ;
            xAy -= Ay [j] ;

        }
        else if (sum_y > 1) {
            j = ind_bad_y [num_bad_y - 1] ;
            if(!IN_RANGE(j, 0, n))
            {
                printf("\n\nBADBADBAD\n\n") ;
                fflush(stdout) ;
            }

            y[j] = 0 ;
            num_bad_y -= 1 ;
            sum_y-- ;

            Az_pm_col_of_A(Ay, j, MINUS, problem) ;
            xAy -= Ax [j] ;
        }
        else
        {
            /* Terrible situation. Just find two non-neighboring vertices and
               start there */
            terribleSituation = true ;
            // printf("\nEntered terrible situation") ;
            // printf("\nGraph name: %s \n", problem->name) ;
            // fflush(stdout) ;
            VSPLONG *Ap = problem->Graph->Ap ;
            VSPINT *Ai = problem->Graph->Ai ;
            VSPINT p ;
            VSPLONG k ;
            bool complete_graph = true ;

            /* Iterate through all the vertices */
            for(i = 0 ; i < n ; i++ )
            {
                /* Vertex i is connected to everything if it has n - 1 ones */
                if(Ap[i+1] - Ap[i] < n - 1)
                {
                    complete_graph = false ;

                    /* Create a hash map of neighbors of i */
                    VSPINT *neighbors = calloc(n, sizeof(*neighbors)) ;
                    
                    VSPLONG end = Ap[i + 1] ;
                    for(k = Ap[i] ; k < end ; k++) neighbors[ Ai[k] ] = 1 ;

                    for(k = 0 ; k < n ; k++)
                    {
                        if(neighbors[k] == 0 && k != i) break ;
                    }
                    #ifdef DEBUG
                        if(k >= n) Error("Problem", DBG_INFO) ;
                    #endif

                    /* Put i in x and k in y */
                    for(p = 0 ; p < n ; p++)
                    {
                        x[p] = 0 ;
                        y[p] = 0 ;
                    }
                    x[i] = 1 ;
                    y[k] = 1 ;
                    X->weight = c[i] ;
                    Y->weight = c[k] ;
                    Compute_Az(X->point, X->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                    Compute_Az(Y->point, Y->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                    printf("\nCompute Az called!\n") ;
                    fflush(stdout) ;

                    free(neighbors) ;

                    break ;
                }
            }
            if(complete_graph)
            {
                printf("\n\nERROR: Graph is complete.") ;
                printf("\nGraph name: %s\n\n", problem->name) ;
                exit(1) ;
            }
            break ;
        }
    }



    /* TODO: (this is a shortcut. Just compute weight above) */
    /* Compute final weight */
    VSPINT weight_x = 0 ;
    VSPINT weight_y = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        weight_x += x [i] * c [i] ;
        weight_y += y [i] * c [i] ;
    }
    X->weight = weight_x ;
    Y->weight = weight_y ;

    /* Free memory */
    free(ind_bad_x) ;
    free(ind_bad_y) ;


    #ifdef DEBUG
        check_states(X, Y, problem) ;

        /* Check that xAy = 0 */
        VSPINT xAy_dbg = 0, yAx_dbg = 0 ;
        for (i = 0 ; i < n ; i++) {
            xAy_dbg += x [i] * Ay [i] ;
            yAx_dbg += y [i] * Ax [i] ;
        }
        if ( xAy_dbg != 0 || yAx_dbg != 0 ) {
            Error("Contraint not satisfied", DBG_INFO) ;
        }

        /* Check that the objective has improved */
        VSPINT obj_new = 0 ;

        for (i = 0 ; i < n ; i++) {
            obj_new += c [i] * ( x [i] + y [i] ) - gamma * x [i] * Ay [i] ;
        }

        if (obj_new < obj_old && terribleSituation == false)
        {
            Error("Objective worsened when finding separator", DBG_INFO) ;
        }
    #endif

    improve(X, Y, problem) ;

    #ifdef DEBUG
        check_states(X, Y, problem) ;

        /* Check that xAy = 0 */
        xAy_dbg = 0 ;
        yAx_dbg = 0 ;
        for (i = 0 ; i < n ; i++) {
            xAy_dbg += x [i] * Ay [i] ;
            yAx_dbg += y [i] * Ax [i] ;
        }
        if ( xAy_dbg != 0 || yAx_dbg != 0 ) {
            Error("Contraint not satisfied", DBG_INFO) ;
        }

        /* Check that the objective has improved */
        obj_new = 0 ;

        for (i = 0 ; i < n ; i++) {
            obj_new += c [i] * ( x [i] + y [i] ) - gamma * x [i] * Ay [i] ;
        }

        if (obj_new < obj_old && terribleSituation == false) {
            Error("Objective worsened when improving point", DBG_INFO) ;
        }

        VSPFLOAT objout = compute_obj(X, Y, problem, problem->Graph->max_vtx_cost) ;
        if(objout < objin && terribleSituation == false)
        {
            Error("Objective worsened after improve", DBG_INFO) ;
        }
    #endif

    TIMER_END ;
    PERFORMANCE_UPDATE_TIME(Relaxed_To_VSP) ;

    return ;
}

/* Assumes (x, y) is fully binary */
void improve
(
    State *X,
    State *Y,
    Problem *problem
)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;
    #endif


    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;
    VSPINT *Ax = X->Az->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT *c = problem->Graph->c ;
    VSPINT u = problem->u ;


    #ifdef DEBUG
        check_states(X, Y, problem) ;
        VSPINT xAy1 = 0 ;
        VSPINT xAy2 = 0 ;
        VSPINT k ;
    #endif

    /* Make a pass through all nodes in the separator and see if any can
       be brought into x or y */
    VSPINT num_improves = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        VSPINT ci = c[i] ;
        if(x[i] == 0 && y[i] == 0)
        {
            /* If i is not connected to either side */
            if(Ax[i] == 0 && Ay[i] == 0)
            {
                #ifdef DEBUG
                        check_states(X, Y, problem) ;
                #endif
                /* i has no connections to x or y. Bring it into the set with
                   smaller weight */
                /* TODO long : Check Jamse's code. When u_a and u_b are
                   different, we want to to things a bit differently here */
                if( Y->weight + ci <= u && Y->weight < X->weight )
                {
                    /* Move i to y */
                    y[i] = 1 ;
                    Y->weight += ci ;
                    // Compute_Az(Y->point, Y->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                    Az_pm_col_of_A(Y->Az->int_part, i, PLUS, problem) ;
                    // printf("\nImprove worked 1\n") ;
                    #ifdef DEBUG
                        check_states(X, Y, problem) ;
                        xAy1 = 0 ;
                        xAy2 = 0 ;
                        for(k = 0 ; k < n ; k++)
                        {
                            xAy1 += x[k] * Ay[k] ;
                            xAy2 += y[k] * Ax[k] ;
                        }
                        if(xAy1 != 0 || xAy1 != xAy2)
                        {
                            Error("Not a separator", DBG_INFO) ;
                        }
                    #endif
                    num_improves++ ;
                }
                else if( X->weight + ci <= u)
                {
                    /* Move i to x */
                    x[i] = 1 ;
                    X->weight += ci ;
                    // Compute_Az(X->point, X->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                    Az_pm_col_of_A(X->Az->int_part, i, PLUS, problem) ;
                    // printf("\nImprove worked 2\n") ;
                    #ifdef DEBUG
                        check_state(X, problem) ;
                        check_state(Y, problem) ;
                        xAy1 = 0 ;
                        xAy2 = 0 ;
                        for(k = 0 ; k < n ; k++)
                        {
                            xAy1 += x[k] * Ay[k] ;
                            xAy2 += y[k] * Ax[k] ;
                        }
                        if(xAy1 != 0 || xAy1 != xAy2)
                        {
                            Error("Not a separator", DBG_INFO) ;
                        }
                    #endif
                    num_improves++ ;
                }
            }
            /* If i is connected to A, and it fits */
            else if(Ay[i] == 0 && X->weight + ci <= u)
            {
                /* Move i to x */
                x[i] = 1 ;
                X->weight += ci ;
                // Compute_Az(X->point, X->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                Az_pm_col_of_A(X->Az->int_part, i, PLUS, problem) ;
                // printf("\nImprove worked 3\n") ;
                #ifdef DEBUG
                    check_state(X, problem) ;
                    check_state(Y, problem) ;
                        xAy1 = 0 ;
                        xAy2 = 0 ;
                        for(k = 0 ; k < n ; k++)
                        {
                            xAy1 += x[k] * Ay[k] ;
                            xAy2 += y[k] * Ax[k] ;
                        }
                        if(xAy1 != 0 || xAy1 != xAy2)
                        {
                            Error("Not a separator", DBG_INFO) ;
                        }
                #endif
                num_improves++ ;
            }
            /* If i is connected to B, and it fits */
            else if(Ax[i] == 0 && Y->weight + ci <= u)
            {
                /* Move i to y */
                y[i] = 1 ;
                Y->weight += ci ;
                // Compute_Az(Y->point, Y->Az, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->n) ;
                Az_pm_col_of_A(Y->Az->int_part, i, PLUS, problem) ;
                // printf("\nImprove worked 4\n") ;
                #ifdef DEBUG
                    check_state(X, problem) ;
                    check_state(Y, problem) ;
                        xAy1 = 0 ;
                        xAy2 = 0 ;
                        for(k = 0 ; k < n ; k++)
                        {
                            xAy1 += x[k] * Ay[k] ;
                            xAy2 += y[k] * Ax[k] ;
                        }
                        if(xAy1 != 0 || xAy1 != xAy2)
                        {
                            Error("Not a separator", DBG_INFO) ;
                        }
                #endif
                num_improves++ ;
            }
        }
    }
    #ifdef DEBUG
        check_states(X, Y, problem) ;
        xAy1 = 0 ;
        xAy2 = 0 ;
        for(k = 0 ; k < n ; k++)
        {
            xAy1 += x[k] * Ay[k] ;
            xAy2 += y[k] * Ax[k] ;
        }
        if(xAy1 != 0 || xAy1 != xAy2)
        {
            Error("Not a separator", DBG_INFO) ;
        }
    #endif
}
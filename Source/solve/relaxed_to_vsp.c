#include "core.h"
#include "solve_internal.h"
#include "utils.h"
#include <stdio.h>
#include "profiler.h"

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
    VSPINT n = problem->graph->n ;
    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;
    VSPINT *Ax = X->Az->int_part ;
    VSPINT *Ay = Y->Az->int_part ;
    VSPINT *c = problem->graph->c ;
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
                    // Compute_Az(Y->point, Y->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
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
                    // Compute_Az(X->point, X->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
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
                // Compute_Az(X->point, X->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
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
                // Compute_Az(Y->point, Y->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
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
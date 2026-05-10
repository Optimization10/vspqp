#include "VSPQP.h"
#include "core.h"
#include "solve_internal.h"
#include "mc.h"
#include "utils.h"
#include "VS.h"

static void set_manual_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
)
{
    VSPINT n = problem->graph->n ;
    VSPINT i ;

    for(i = 0 ; i < n ; i++) {
        X->point->int_part[i] = 0 ;
        Y->point->int_part[i] = 0 ;
    }
    X->point->int_part[0] = 1 ;
    Y->point->int_part[0] = 1 ;

    X->point->ind_frac = EMPTY ;
    X->point->val_frac = 0 ;

    Y->point->ind_frac = EMPTY ;
    Y->point->val_frac = 0 ;

    compute_weight(X->point, &(X->weight), problem->graph->c) ;
    compute_weight(Y->point, &(Y->weight), problem->graph->c) ;

    Compute_Az(X->point, X->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
    Compute_Az(Y->point, Y->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;

    #ifdef DEBUG
        VSPFLOAT obj_dbg ;
        obj_dbg = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj_dbg ;
    #endif

    return ;
}



static void set_random_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
)
{
    VSPINT i, k, temp, xwgt, ywgt ;
    VSPINT n = problem->graph->n ;

    VSPINT *x = X->point->int_part ;
    VSPINT *y = Y->point->int_part ;

    /* Set x and y to 0 */
    for(i = 0 ; i < n ; i++)
    {
        x[i] = 0 ;
        y[i] = 0 ;
    }

    /* Get a random permutation of vertices */
    VSPINT *perm_1 = malloc(n * sizeof(*perm_1)) ;
    VSPINT *perm_2 = malloc(n * sizeof(*perm_2)) ;
    for(i = 0 ; i < n ; i++)
    {
        perm_1[i] = i ;
        perm_2[i] = i ;
    }


    for(i = 0 ; i < n ; i++)
    {
        /* Find a random k in [0, n-1] */
        k = rand() % n ;

        /* Exchange perm_1[i]  with perm_1[k] */
        temp = perm_1[i] ;
        perm_1[i] = perm_1[k] ;
        perm_1[k] = temp ;

        /* Find a random k in [0, n-1] */
        k = rand() % n ;

        /* Exchange perm_2[i]  with perm_2[k] */
        temp = perm_2[i] ;
        perm_2[i] = perm_2[k] ;
        perm_2[k] = temp ;

    }

    /* Fill out x and y randomly */
    VSPINT target_weight = ceil(problem->graph->total_vtx_cost * options->random_point_weight_proportion) ;
    target_weight = MIN(target_weight, problem->u) ;
    target_weight = MAX(target_weight, 1) ;

    xwgt = 0 ;
    ywgt = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        k = perm_1[i] ;
        if(xwgt + problem->graph->c[k] <= target_weight)
        {
            x[k] = 1 ;
            xwgt += problem->graph->c[k] ;
        }
        else
        {
            x[k] = 0 ;
        }


        k = perm_2[i] ; 
        if(ywgt + problem->graph->c[k] <= target_weight)
        {
            y[k] = 1 ;
            ywgt += problem->graph->c[k] ;
        }
        else
        {
            y[k] = 0 ;
        }
    }

    /* Update everything */
    X->weight = xwgt ;
    X->point->ind_frac = EMPTY ;
    X->point->val_frac = 0 ;
    Compute_Az(X->point, X->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;

    Y->weight = ywgt ;
    Y->point->ind_frac = EMPTY ;
    Y->point->val_frac = 0 ;
    Compute_Az(Y->point, Y->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;

    #ifdef DEBUG
        check_states(X, Y, problem) ;
    #endif

    /* Free memory */
    free(perm_1) ;
    free(perm_2) ;

}


/*  Start at the centroid and do an iteration of mountain climbing.
    The reason for mc here is that the mc only accepts, as input, almost binary
    points.

    In principle, the function does the following:
        Let the centroid be (x_c, y_c). Do:
        - Obtain x_star by maximizing f(x  , y_c) over x
        - Obtain y_star by maximizing f(x_c, y  ) over y
        - If f(x_star, y_star) is greater than both
          f(x_star, y) and f(x, y_star), then return x_star, y_star
        - Else, obtain x_star2 by maximizing f(x, y_star) over x, and return
          (x_star2, y_star)

    The centroid is (x0, y0) with
        x0 = f * 1,    y0 = f * 1,
    where 1 is the vector of all 1's, and f = bal / sum_i c_i.
    At the centroid, the gradient of f with respect to x is
    c - gamma A * y0 = c - gamma * f * A * 1
    The ith entry of A * 1 can be computed efficiently as Ap[i+1] - Ap[i]
     */
static void set_centroid
(
    State             *X,
    State             *Y,
    Problem     *problem,
    Options     *options
)
{
    VSPINT i ;
    VSPINT A1i ;
    VSPINT n = problem->graph->n ;
    VSPINT u = problem->u ;
    VSPFLOAT f ;
    VSPINT gamma = problem->graph->max_vtx_cost ;
    VSPLONG *Ap = problem->graph->Ap ;
    VSPINT *c = problem->graph->c ;

    /* Allocate memory */
    VSPFLOAT   *d = calloc(n, sizeof(*d)) ;
    VSPFLOAT *y_f = calloc(n, sizeof(*y_f)) ;

    /* Initialize the centroid */
    f = options->centroidVal ;

    /* Compute d, the gradient at the centroid */
    for(i = 0 ; i < n ; i++) {
        A1i = Ap[i+1] - Ap[i] + 1 ;       /* The 1 is for the diagonal */
        d[i] = c[i] - gamma * f * A1i ;
    }

    #ifdef DEBUG
        /* Compute objective */
        VSPFLOAT obj_dbg = 0 ;
        for(i = 0 ; i < n ; i++) {
            obj_dbg += f * d[i] + f * c[i] ;
        }
    #endif

    /* Optimize over y and compute the objective */
    VSPFLOAT obj_y = LP_solver_0(y_f, d, c, 1, u, n, options) ;
    obj_y += problem->u ;

    /* Update the state of y */
    array_to_point(y_f, Y->point) ;
    Compute_Az(Y->point, Y->Az, problem->graph->Ai, problem->graph->Ap, problem->graph->n) ;
    compute_weight(Y->point, &Y->weight, c) ;

    /* Check objective at diagonal step. Diagonal step has the form (y, y) */
    VSPFLOAT obj_diag = compute_obj(Y, Y, problem, problem->graph->max_vtx_cost) ;
    if(obj_diag > obj_y)
    {
        /* Take the diagonal step by making X identical to Y */
        copy_state(Y, X) ;
    }
    else
    {
        /* Optimize over X, to move it away from the centroid */
        LP_solver_4(X, Y->Az, c, u, problem->graph->max_vtx_cost, problem,options, 'x');
    }


    #ifdef DEBUG
        /* Compute objective */
        obj_dbg = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj_dbg ;

        /* Check output */
        check_states(X, Y, problem) ;
    #endif


    /* Free memory */
    free(y_f) ;
    free(d) ;


    return ;
}

// If separator from part array is empty, and the shores are feasible
// at the finest level, then only the int parts of X and Y
// are set. Ax and Ay are not set
void Part_To_States
(
    int        *part,
    State      *X,
    State      *Y,
    Problem    *problem,
    bool        empty_sep
)
{
    VSPINT   i ;
    VSPINT   side ;
    VSPINT   u_finest ;
    VSPINT   wgts_part[3] = {0, 0, 0} ;
    VSPINT  *x  = X->point->int_part ;
    VSPINT  *y  = Y->point->int_part ;
    VSPINT   n  = problem->graph->n ;
    VSPINT  *c  = problem->graph->c ;
    VSPINT   u  = problem->u ;
    VSPINT  *Ai = problem->graph->Ai ;
    VSPLONG *Ap = problem->graph->Ap ;

    // ====================================================================== //
    // Populate states according to partition
    // ====================================================================== //
    // Set the points and weights
    X->weight = 0 ;
    Y->weight = 0 ;
    for(i = 0 ; i < n ; i++)
    {
        side = part[i] ;
        wgts_part[side] += c[i] ;

        if(side == 0)
        {
            y[i] = 0 ;
            if(X->weight + c[i] <= u)
            {
                x[i] = 1 ;
                X->weight += c[i] ;
            }
            else x[i] = 0 ;
        }
        else if(side == 1)
        {
            x[i] = 0 ;
            if(Y->weight + c[i] <= u)
            {
                y[i] = 1 ;
                Y->weight += c[i] ;
            }
            else y[i] = 0 ;
        }
        else
        {
            VSP_ASSERT(side == 2) ; 
            x[i] = 0 ;
            y[i] = 0 ;
        }
    }

    // Set fractional entries of the points to 0
    X->point->ind_frac = EMPTY ;
    Y->point->ind_frac = EMPTY ;
    X->point->val_frac = 0 ;
    Y->point->val_frac = 0 ;

    u_finest = Get_u_Finest(problem) ;
    if(empty_sep && wgts_part[2] == 0 && wgts_part[0] <= u_finest && wgts_part[1]<=u_finest)
    {
        for(i = 0 ; i < n ; i++)
        {
            x[i] = 0 ;
            y[i] = 0 ;
            if(part[i] == 0)      x[i] = 1 ;
            else if(part[i] == 1) y[i] = 1 ;
        }
        X->weight = wgts_part[0] ;
        Y->weight = wgts_part[1] ;
        if(empty_sep)
        {
            printf("inside %s\n", G_Current_Problem_Name) ;
        }
    }
    else
    {
        // Set Ax and Ay
        Compute_Az(X->point, X->Az, Ai, Ap, n) ;
        Compute_Az(Y->point, Y->Az, Ai, Ap, n) ;
    }


#ifdef DEBUG
    check_states(X, Y, problem) ;
#endif

}

void SetUpperBound
(
    Problem *problem,
    Options *opts,
    int quadratic
) ;

static void set_bisection_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
)
{
    SET_FUNCTION_COUNTER ;

    VSPINT   i ;
    VSPINT   u    = problem->u ;
    VSPINT   n    = problem->graph->n ;
    VSPINT  *c    = problem->graph->c ;
    VSPINT  *Ai   = problem->graph->Ai ;
    VSPINT  *Ae   = problem->graph->Ae ;
    VSPINT  *part = malloc(n * sizeof(*part)) ;
    VSPINT  *x    = X->point->int_part ;
    VSPINT  *y    = Y->point->int_part ;
    VSPLONG *Ap   = problem->graph->Ap ;

    // Set parameters for the bisection
    VSparm ParmStruc, *Parm ;
    Parm = &ParmStruc ;
    vs_default (Parm) ;
    Parm->initial_project = options->initial_project ;

    // Get partition
    VSPINT m = (problem->graph->total_vtx_cost) / 2 ;
    vs(part, Ai, Ap, Ae, c, n, m,Parm) ;


    if(options->ub_based_on_bisection)
    {
        Problem *P = problem ;
        VSPFLOAT ratioUBLimit ;
        VSPFLOAT ratioUBLimit_Old ;
        VSPINT   wgts[3] = {0, 0, 0} ;

        // Calculate new upper bound at coarsest level
        for(i = 0 ; i < n ; i++) wgts[ part[i] ] += c[i] ;
        ratioUBLimit = MAX(wgts[0], wgts[1]) / (VSPFLOAT) (P->graph->total_vtx_cost) ;
        ratioUBLimit = MAX(ratioUBLimit, options->ratioUBLimit) ;

        // printf("Coarsest ub: %f\n", ratioUBLimit) ;

        ratioUBLimit_Old = options->ratioUBLimit ;
        options->ratioUBLimit = ratioUBLimit ;
        while(P != NULL)
        {
            SetUpperBound(P, options, 1) ;
            P = P->parent_problem ;
        }

        // Restore
        options->ratioUBLimit = ratioUBLimit_Old ;
    }

    // Convert partition to states
    Part_To_States(part, X, Y, problem, true) ;

    // Cleanup
    free(part) ;
}

// May return a point that's infeasible for the coarsest level, but
// has zero sep and is feasible for the finest level
void set_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
)
{
    switch(options->initial_guess)
    {
        case centroid_initial_guess :
            set_centroid(X, Y, problem, options) ;
            break ;
        case manual_initial_guess :
            set_manual_initial_guess(X, Y, problem, options) ;
            break ;
        case random_initial_guess :
            set_random_initial_guess(X, Y, problem, options) ;
            break ;
        case bisection_initial_guess :
            set_bisection_initial_guess(X, Y, problem, options) ;
            break ;
    }
    #ifdef DEBUG
        VSPFLOAT obj = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj ;
    #endif
}
#include "solve_internal.h"
#include "mc_api.h"


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
    }
    #ifdef DEBUG
        VSPFLOAT obj = compute_obj(X, Y, problem, problem->graph->max_vtx_cost) ;
        (void) obj ;
    #endif

    return ;
}
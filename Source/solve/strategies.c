#include "VSPQP.h"
#include "solve_internal.h"
#include "utils.h"
#include "greedy_flips.h"
#include "FM.h"
#include "VS.h"
#include "mc.h"

// Run a strategy on a cloned point. If it beats the original, swap.
// Returns true if the clone was better.
// ratio is the ratio between the two objectives. ratio can be NULL
bool Try_And_Keep_If_Better
(
    State    *X,
    State    *Y,
    Problem  *problem,
    void    (*strategy)(State*, State*, Problem*, Options*),
    Options  *options,
    VSPFLOAT *ratio    
)
{
    VSPINT total_vtx_cost ;
    VSPFLOAT sep_new, sep_orig ;
    State *Xc = Clone_State(X);
    State *Yc = Clone_State(Y);
    VSPINT gamma = problem->graph->max_vtx_cost;

    strategy(Xc, Yc, problem, options);

    VSPFLOAT obj_orig = compute_obj(X,  Y,  problem, gamma);
    VSPFLOAT obj_new  = compute_obj(Xc, Yc, problem, gamma);

    bool better = (obj_new >= obj_orig);
    if (better)
    {
        Swap_Mem(X, Xc, sizeof(State));
        Swap_Mem(Y, Yc, sizeof(State));
    }

    Destroy_State(Xc);
    Destroy_State(Yc);

    // Set ratio
    total_vtx_cost = problem->graph->total_vtx_cost ;

    sep_orig = total_vtx_cost - obj_orig ;
    if(sep_orig == 0) sep_orig = 1 ;

    sep_new = total_vtx_cost - obj_new ;
    if(sep_new == 0) sep_new = 1 ;

    if(ratio != NULL) *ratio = sep_new / sep_orig ;

    return better;
}

void strategy_centroid_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = centroid_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}

void strategy_bisection_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = bisection_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    if(X->weight + Y->weight == problem->graph->total_vtx_cost)
    {
        printf("intermediate-empty-sep!!!!\n") ;
    }
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}

void strategy_random_restart
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    initial_guess tmp = options->initial_guess ;
    options->initial_guess = random_initial_guess ;
    set_initial_guess(X, Y, problem, options);
    Optimize(X, Y, problem, options);
    Round_Down_Fractions(X, Y, problem);
    options->initial_guess = tmp ;
}


void strategy_explore_bisection
(
    State *X, State *Y, Problem *problem, Options *options
)
{
    // Convert the current separator to an x array for vs()
    VSPINT   i ;
    VSPINT  *x  = X->point->int_part ;
    VSPINT  *y  = Y->point->int_part ;
    VSPINT   n  = problem->graph->n ;
    VSPINT   m  = (problem->graph->total_vtx_cost) / 2 ;
    VSPINT  *c  = problem->graph->c ;
    VSPINT  *Ai = problem->graph->Ai ;
    VSPINT  *Ae = problem->graph->Ae ;
    VSPLONG *Ap = problem->graph->Ap ;
    VSPINT   idx_frac_x = X->point->ind_frac ;
    VSPINT   idx_frac_y = Y->point->ind_frac ;
    VSPFLOAT gamma      = problem->graph->max_vtx_cost ;

    int *part = malloc(n * sizeof(*part)) ;

    for(i = 0 ; i < n ; i++)
    {
        if     (x[i]) part[i] = 0 ;
        else if(y[i]) part[i] = 1 ;
        else          part[i] = 2 ;
    }

    if(idx_frac_x != EMPTY) part[idx_frac_x] = 0 ;
    if(idx_frac_y != EMPTY) part[idx_frac_y] = 1 ;

    // Call vs
    VSparm ParmStruc, *Parm ;
    Parm = &ParmStruc ;
    vs_default (Parm) ;
    Parm->ngp = 5 ;
    Parm->compute_eig = FALSE ;
    vs(part, Ai, Ap, Ae, c, n, m,Parm) ;

    // Convert vs point to states
    Part_To_States(part, X, Y, problem, false) ;
    free(part) ;

    // Do mc
    Optimize(X, Y, problem, options) ;
}
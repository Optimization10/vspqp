#include "core.h"
#include "solve_internal.h"
#include "mc.h"
#include "FM.h"
#include "greedy_flips.h"
#include "utils.h"
#include "vs.h"




void Best_Of_MC_And_Greedy_Flips
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options,
    VSPINT       numFlipsRuns
)
{
    VSPINT gamma = problem->graph->max_vtx_cost;
    VSPINT *c = problem->graph->c;

    // Clone for flips path
    State *Xf = Clone_State(X);
    State *Yf = Clone_State(Y);

    VSPINT wgtx = State_Get_Wgt_Int_Part(Xf, c);
    VSPINT wgty = State_Get_Wgt_Int_Part(Yf, c);

    // Flips path
    if (wgtx != 0 && wgty != 0)
    {
        Multi_Start_Flips(Xf, Yf, problem, options, numFlipsRuns, 9) ;
    }

    // MC path (in-place on X, Y)
    mountain_climbing(X, Y, problem, options, gamma, true);


    // Compare
    VSPFLOAT objMC, objFlips ;
    objMC    = compute_obj(X, Y   , problem, (VSPFLOAT) gamma) ;
    objFlips = compute_obj(Xf, Yf , problem, (VSPFLOAT) gamma) ;

    bool flipsBetter    = (objFlips > objMC);
    bool flipsSame      = (objFlips == objMC);
    bool flipsBetterBal = (IMBALANCE(Xf->weight, Yf->weight) < IMBALANCE(X->weight, Y->weight));

    if (wgtx != 0 && wgty != 0 && (flipsBetter || (flipsSame && flipsBetterBal)))
    {
        Swap_Mem(X, Xf, sizeof(State));
        Swap_Mem(Y, Yf, sizeof(State));
    }

    Destroy_State(Xf);
    Destroy_State(Yf);
}

static void Determine_Num_GR
(
    VSPINT             n,
    Options *options,
    VSPINT        *numGR  // output
)
{
    *numGR = 0 ;
    if(n <= options->max_num_vert_small_graph)
    {
        *numGR = options->num_gr_for_small_graphs ;
    }
    else
    {
        *numGR = options->num_gr_for_large_graphs ;
    }
    return ;
}


// Single an early exit if f(x, y) >= c(V) - alpha
// The logic is that if (X, Y) represent a vertex separator
// partition, and the separator only contains a node of minimal weight, then
// X, Y is globally optimal, or nearly globally optimal
// Note that we are not taking into account the penalty term
static bool Early_Exit(State *X, State *Y, Problem *P)
{
    VSPFLOAT obj = compute_obj(X, Y, P, P->graph->max_vtx_cost) ;
    bool ret = false ;
    if(lteq(P->graph->total_vtx_cost - P->graph->min_vtx_cost, obj))
    {
        ret = true ;
    }
    return ret ;
}


static void Gamma_Reduction_Loop
(
    State       *X,
    State       *Y,
    Problem     *problem,
    Options     *options
)
{
    VSPINT i, k = 0, num_gr;
    VSPINT gamma_org = problem->graph->max_vtx_cost ;
    VSPFLOAT gamma, f_x_y, f_xstar_ystar;

    Determine_Num_GR(problem->graph->n, options, &num_gr);

    for (i = 0; i < num_gr; i++)
    {
        if (k == 0)
        {
            gamma = minreduce(X, Y, options, problem) - 1E-3;
            if (gamma < 0) break;
        }
        else gamma = gamma / 1.3;

        State *Xc = Clone_State(X);
        State *Yc = Clone_State(Y);

        mountain_climbing(Xc, Yc, problem, options, gamma, false);
        mountain_climbing(Xc, Yc, problem, options, gamma_org, true);
        Round_Down_Fractions(Xc, Yc, problem) ;
        Greedy_Flips(Xc, Yc, problem, 9, (problem->level <= 2), false);

        VSPINT np;
        if (problem->level >= problem->coarsest_level_num / 2) np = 3;
        else np = 2;
        FM_Metis(Xc, Yc, problem, options, np);

        f_x_y         = compute_obj(X,  Y,  problem, gamma_org);
        f_xstar_ystar = compute_obj(Xc, Yc, problem, gamma_org);

        if (f_xstar_ystar > f_x_y + fabs(f_x_y) * 1E-13)
        {
            Swap_Mem(X, Xc, sizeof(State));
            Swap_Mem(Y, Yc, sizeof(State));
            k = 0;

            if (Early_Exit(X, Y, problem))
            {
                Destroy_State(Xc);
                Destroy_State(Yc);
                return;
            }
        }
        else k++;

        Destroy_State(Xc);
        Destroy_State(Yc);
    }
}


void Optimize
(
    State   *X,
    State   *Y,
    Problem *problem,
    Options *options
)
{
    SET_FUNCTION_COUNTER ;
#ifdef DEBUG
        check_states(X, Y, problem);
#endif
    VSPFLOAT *x_flt, *y_flt, *gx, *gy, *c_flt, xcost, ycost ;
    VSPINT    n_passes;
    VSPINT    n     = problem->graph->n ;
    VSPINT   *c     = problem->graph->c ;
    VSPINT    l     = problem->level ;
    VSPINT    L     = problem->coarsest_level_num ;
    VSPINT    gamma = problem->graph->max_vtx_cost ;
    VSPINT   *Ai    = problem->graph->Ai ;
    VSPLONG  *Ap    = problem->graph->Ap ;
    VSPINT    u     = problem->u ;

    /* ====================================================================== */
    // Phase 0: Gradient Projection
    /* ====================================================================== */
    if(options->vsp_grad_proj)
    {
        /* Step 1: Convert x, y, c to floats, and compute gx and gy */
        x_flt = Point_To_Array(X->point) ;
        gx    = Az_To_Grad(c, gamma, Y->Az) ;

        y_flt = Point_To_Array(Y->point) ;    
        gy    = Az_To_Grad(c, gamma, X->Az) ;

        c_flt = Array_Int_To_Flt(c, n) ;

        /* Step 2: Do gradient projection */
        VSP_Grad_Proj
        (
            x_flt, gx, y_flt, gy, c_flt, Ai, Ap, u, n, gamma, options->vsp_ngp
        ) ;

        /* Step 3: Move x_flt, y_flt to almost binary points */
        int print_level = 0 ;

        VSPINT j1, j2 ;
        j1 = VSPboundary(x_flt, c_flt, gx, gy, n, Ai, Ap, gamma, print_level) ;
        j2 = VSPboundary(y_flt, c_flt, gy, gx, n, Ai, Ap, gamma, print_level) ;

        /* Step 4: Update the X and Y states with x_flt and y_flt */
        xcost = Move_To_Vertex(x_flt, y_flt, gx, gy, Ai, Ap, c_flt, u, gamma, n, &j1) ;
        ycost = Move_To_Vertex(y_flt, x_flt, gy, gx, Ai, Ap, c_flt, u, gamma, n, &j2) ;

        Array_To_State(x_flt, xcost, gy, X, problem, j1) ;
        Array_To_State(y_flt, ycost, gx, Y, problem, j2) ;

        /* Step 5: Cleanup */
        free(x_flt) ;
        free(y_flt) ;
        free(gx) ;
        free(gy) ;
        free(c_flt) ;
    }

    /* ====================================================================== */
    // Phase 1: Mountain climbing
    /* ====================================================================== */
    mountain_climbing(X, Y, problem, options, gamma, true) ;
    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 2: Greedy flips
    Round_Down_Fractions(X, Y, problem) ;
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2), false);

    // Phase 3: FM
    n_passes = (l >= L / 2) ? 3 : 2 ;
    FM_Metis(X, Y, problem, options, n_passes);
    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    // Phase 4: MC
    mountain_climbing(X, Y, problem, options, gamma, true) ;

    // Phase 5: Gamma reductions
    Gamma_Reduction_Loop(X, Y, problem, options);

    // Phase 6: Final cleanup flips
    Round_Down_Fractions(X, Y, problem) ;
    Greedy_Flips(X, Y, problem, 5, true, true);

    CLEANUP:

    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    return;
}
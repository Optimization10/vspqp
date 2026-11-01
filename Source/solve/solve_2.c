#include "VSPQP.h"
#include "core.h"
#include "solve_internal.h"
#include "coarsen.h"
#include "utils.h"
#include "profiler.h"
#include "greedy_flips.h"
#include "FM.h"
#include <stdint.h>
#include <string.h>

/* ========================================================================== */
// VSP_Solve_2: METIS-style branched coarsening
/* ========================================================================== */
// The finest graph is coarsened branch_shared_levels times, giving the root
// of the branches. The root is then coarsened num_branches times
// independently, each time down to fewer than branch_node_cutoff vertices.
// Each branch is solved at its coarsest level and refined back up to the
// root. The best root solution is refined up to the finest level.
//
// Branch 0 visits the vertices in index order when matching, like VSP_Solve,
// so with one branch (and branch_full_solves >= 1) VSP_Solve_2 gives
// VSP_Solve's result for node_cutoff = branch_node_cutoff. The other branches
// visit the vertices in a new uniformly random order at each level, then
// renumber the coarse vertices so that each coarse graph keeps the vertex
// order of its finer graph.
//
// The first branch_full_solves branches get the full Solve_Level on the levels
// they solve on their own (their coarsest level up to and including the root).
// The others get a cheap solve there: greedy flips and FM only. If the winning
// root solution comes from a cheap branch, it gets the full Solve_Level at the
// root first. Every finer level always gets the full Solve_Level.


/* -------------------------------------------------------------------------- */
// Random visit orders
/* -------------------------------------------------------------------------- */

// splitmix64. Kept apart from rand() so the visit orders don't depend on the
// rand() calls made by the solvers, and don't change them either.
static uint64_t Rng_Next(uint64_t *state)
{
    uint64_t z = (*state += 0x9E3779B97F4A7C15ull) ;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull ;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull ;
    return z ^ (z >> 31) ;
}

// perm = a uniformly random permutation of 0, 1, ..., n-1
static void Random_Order(VSPINT *perm, VSPINT n, uint64_t *rng)
{
    for(VSPINT i = 0 ; i < n ; i++) perm[i] = i ;
    for(VSPINT i = n - 1 ; i > 0 ; i--)
    {
        VSPINT j = (VSPINT) (Rng_Next(rng) % (uint64_t) (i + 1)) ;
        SWAP(perm[i], perm[j], VSPINT) ;
    }
}

// Matching numbers the coarse vertices in the order their matches are made,
// so after a random visit order the coarse graph's numbering is random too.
// Renumber the coarse vertices by their lowest-numbered fine vertex, so the
// coarse graph keeps the fine graph's vertex order (METIS does the same).
// newid is workspace of size at least M->cn.
static void Renumber_Coarse_Vertices
(
    Matching *M,
    VSPINT    n,     // number of fine vertices
    VSPINT   *newid
)
{
    VSPINT cv, fv, next = 0 ;
    VSPINT *matchmap    = M->matchmap ;
    VSPINT *invmatchmap = M->invmatchmap ;

    for(cv = 0 ; cv < M->cn ; cv++) newid[cv] = EMPTY ;

    for(fv = 0 ; fv < n ; fv++)
    {
        cv = matchmap[fv] ;
        if(newid[cv] == EMPTY)
        {
            newid[cv]         = next ;
            invmatchmap[next] = fv ;
            next++ ;
        }
        matchmap[fv] = newid[cv] ;
    }

    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count++ ;

        if(next != M->cn) ERR("Wrong number of coarse vertices") ;
        for(cv = 0 ; cv < M->cn ; cv++)
        {
            if(matchmap[ invmatchmap[cv] ] != cv) ERR("Error in invmatchmap") ;
        }
        for(fv = 0 ; fv < n ; fv++)
        {
            // Every vertex maps where its match does
            if(matchmap[ M->match[fv] - 1 ] != matchmap[fv])
            {
                ERR("Error in matchmap") ;
            }
        }
    #endif
}

// Clear a matching so the same graph can be matched again
static void Reset_Matching(Matching *M, VSPINT n)
{
    memset(M->match,       0, n * sizeof(*M->match)) ;
    memset(M->matchmap,    0, n * sizeof(*M->matchmap)) ;
    memset(M->invmatchmap, 0, n * sizeof(*M->invmatchmap)) ;
    memset(M->matchtype,   0, n * sizeof(*M->matchtype)) ;
    M->cn = 0 ;
}


/* -------------------------------------------------------------------------- */
// Coarsening
/* -------------------------------------------------------------------------- */

// Coarsen P at most max_levels times (EMPTY: no limit) and return the last
// problem kept. Uses the stopping rules of Problem_Recursive_Coarsen, with
// branch_node_cutoff in place of node_cutoff. *bottom is set if a stopping
// rule ended the coarsening, so the returned problem is the coarsest level.
// If rng is NULL vertices are matched in index order, otherwise in a random
// order: on every level, or only on the first if branch_permute_once is set.
static Problem *Coarsen
(
    Problem   *P,
    Options   *opts,
    VSPINT     max_levels,
    uint64_t  *rng,
    VSPINT     N0,         // n at the finest level
    VSPINT     u0,         // u at the finest level
    VSPINT   **sing_sol,   // Output: singleton short-circuit solution, or NULL
    VSPINT    *n_sings,    // Output
    bool      *bottom      // Output
)
{
    VSPINT   l, n ;
    VSPLONG  nz ;
    Problem *coarse_P ;
    VSPINT  *perm  = NULL ;
    VSPINT  *newid = NULL ;
    VSPINT   lo    = N0 - u0 ;

    *sing_sol = NULL ;
    *n_sings  = 0 ;
    *bottom   = true ;

    // Graphs only get smaller, so size the workspace for P
    if(rng != NULL)
    {
        perm  = malloc(P->graph->n * sizeof(*perm)) ;
        newid = malloc(P->graph->n * sizeof(*newid)) ;
    }

    for(l = 0 ; max_levels == EMPTY || l < max_levels ; l++)
    {
        TIMER_START(T_COARSEN) ;
        // Match P
        bool random = (rng != NULL) && (l == 0 || !opts->branch_permute_once) ;
        if(random) Random_Order(perm, P->graph->n, rng) ;
        Graph_Simple_Match
        (
            P->graph, P->matching, P->level, opts, random ? perm : NULL
        ) ;
        if(random)
        {
            Renumber_Coarse_Vertices(P->matching, P->graph->n, newid) ;
        }

        // Coarsen P
        coarse_P = Problem_Coarsen(P, opts) ;

        if(coarse_P->graph->singletons_cost >= lo)
        {
            *sing_sol = Singletons_Solution(coarse_P, N0, u0, n_sings) ;
            if(*sing_sol != NULL)
            {
                P = coarse_P ;
                goto DONE ;
            }
        }

        // Stopping conditions
        n  = coarse_P->graph->n ;
        nz = coarse_P->graph->nz ;
        if(Graph_Density(n, nz) > opts->density_cutoff
           || n < opts->branch_node_cutoff)
        {
            if(Is_Complete_Graph(n, nz)) Destroy_Problem(coarse_P) ;
            else                         P = coarse_P ;
            goto DONE ;
        }
        else if(nz == 0)
        {
            Destroy_Problem(coarse_P) ;
            goto DONE ;
        }

        TIMER_STOP(T_COARSEN) ;
        Timer_Flush(P->level) ;
        P = coarse_P ;
    }
    *bottom = false ;

    DONE:
    P->child_problem = NULL ;
    free(perm) ;
    free(newid) ;
    return P ;
}


/* -------------------------------------------------------------------------- */
// Branched solve
/* -------------------------------------------------------------------------- */

// Cheap level solve for the branches: Phases 2 and 4 of Optimize (greedy
// flips, then FM), without mountain climbing or gamma reduction. FM gives most
// of the full solve's quality, and the parts left out most of its time.
static void Solve_Level_Cheap
(
    State    *X,
    State    *Y,
    Problem  *problem,
    Options  *options
)
{
    #ifdef DEBUG
        check_problem(problem, options) ;
    #endif

    TIMER_START(T_TOTAL) ;
    Round_Down_Fractions(X, Y, problem) ;
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2), false) ;
    Round_Down_Fractions(X, Y, problem) ;
    FM_Metis
    (
        X, Y, problem, options,
        (problem->level >= problem->coarsest_level_num / 2) ? 3 : 2
    ) ;
    TIMER_STOP(T_TOTAL) ;
    Timer_Flush(problem->level) ;
}

static void Solve_Branched_Hierarchy
(
    Problem  *probFinest,
    Options  *options,
    State   **X_0,
    State   **Y_0
)
{
    VSPINT  *sol_sings ;
    VSPINT   n_sings ;
    bool     bottom, global ;
    Problem *root, *coarsest, *P ;
    State   *X, *Y ;

    VSPINT N0       = probFinest->graph->n ;
    VSPINT u0       = probFinest->u ;
    VSPINT shared   = MAX(0, options->branch_shared_levels) ;
    VSPINT branches = MAX(1, options->num_branches) ;

    /*------------------------------------------------------------------------*/
    // Shared levels
    /*------------------------------------------------------------------------*/
    root = Coarsen
    (
        probFinest, options, shared, NULL, N0, u0, &sol_sings, &n_sings, &bottom
    ) ;

    // Reached the coarsest level within the shared levels, so there's nothing
    // to branch from. Solve the single hierarchy as VSP_Solve does.
    if(bottom)
    {
        Set_Upper_Bounds(root, options, options->u_ratio_coarsest, QUADRATIC) ;
        global = Solve_Coarsest
        (
            root, options, sol_sings, n_sings, Solve_Level, X_0, Y_0
        ) ;
        Solve_Up(root, probFinest, global, options, Solve_Level, X_0, Y_0) ;
        free(sol_sings) ;
        return ;
    }

    /*------------------------------------------------------------------------*/
    // Branches
    /*------------------------------------------------------------------------*/
    VSPINT    n_root      = root->graph->n ;
    VSPINT    total       = root->graph->total_vtx_cost ;
    VSPFLOAT  gamma       = root->graph->max_vtx_cost ;
    VSPFLOAT  obj, obj_best = 0 ;
    bool      global_best = false ;
    bool      cheap_best  = false ; // true if the best came from a cheap branch
    VSPINT    L_best      = 0 ;
    VSPINT   *u_best      = malloc((root->level + 1) * sizeof(*u_best)) ;
    State    *X_best      = NULL ;
    State    *Y_best      = NULL ;
    uint64_t  seeds       = (uint64_t) options->branch_seed ;

    for(VSPINT b = 0 ; b < branches ; b++)
    {
        uint64_t     rng          = Rng_Next(&seeds) ;
        bool         cheap        = (b >= options->branch_full_solves) ;
        Level_Solver branch_solve = cheap ? Solve_Level_Cheap : Solve_Level ;

        // Coarsen the root from scratch
        Reset_Matching(root->matching, n_root) ;
        coarsest = Coarsen
        (
            root, options, EMPTY, (b == 0) ? NULL : &rng,
            N0, u0, &sol_sings, &n_sings, &bottom
        ) ;
        Set_Upper_Bounds
        (
            coarsest, options, options->u_ratio_coarsest, QUADRATIC
        ) ;

        // Solve the branch and refine it up to the root
        global = Solve_Coarsest
        (
            coarsest, options, sol_sings, n_sings, branch_solve, &X, &Y
        ) ;
        Solve_Up(coarsest, root, global, options, branch_solve, &X, &Y) ;
        free(sol_sings) ;

        Destroy_All_Problems(root->child_problem) ;
        root->child_problem = NULL ;

        // Keep the larger objective, then the better balance. With an empty
        // separator Az was never computed, but the objective is known.
        obj = global ? total : compute_obj(X, Y, root, gamma) ;
        if(   X_best == NULL || obj > obj_best
           || (obj == obj_best && imbalance(X, Y) < imbalance(X_best, Y_best)))
        {
            if(X_best != NULL)
            {
                Destroy_State(X_best) ;
                Destroy_State(Y_best) ;
            }
            X_best      = X ;
            Y_best      = Y ;
            obj_best    = obj ;
            global_best = global ;
            cheap_best  = cheap ;

            // Coarsening sets the depth, and Set_Upper_Bounds the bounds, of
            // the shared levels for this branch. Save them.
            L_best = root->coarsest_level_num ;
            for(P = root ; P != NULL ; P = P->parent_problem)
            {
                u_best[P->level] = P->u ;
            }
        }
        else
        {
            Destroy_State(X) ;
            Destroy_State(Y) ;
        }

        // An empty separator can't be beaten
        if(lteq(total, obj_best)) break ;
    }

    // Restore the shared levels to how the best branch left them
    for(P = root ; P != NULL ; P = P->parent_problem)
    {
        P->u                  = u_best[P->level] ;
        P->coarsest_level_num = L_best ;
    }
    free(u_best) ;

    /*------------------------------------------------------------------------*/
    // Refine the best root solution up to the finest level
    /*------------------------------------------------------------------------*/
    // A cheap winner hasn't had the full solve at the root yet. (An
    // empty-separator solution has no Az and needs none.)
    if(cheap_best && !global_best)
    {
        Solve_Level(X_best, Y_best, root, options) ;
    }
    Solve_Up
    (
        root, probFinest, global_best, options, Solve_Level, &X_best, &Y_best
    ) ;
    *X_0 = X_best ;
    *Y_0 = Y_best ;
}


solve_result *VSP_Solve_2
(
    VSPINT    n_input,
    VSPLONG  *Ap_input,
    VSPINT   *Ai_input,
    Options  *options,
    char     *name       // Optional. Problem name
)
{
    return VSP_Solve_With
           (
               n_input, Ap_input, Ai_input, options, name,
               Solve_Branched_Hierarchy
           ) ;
}

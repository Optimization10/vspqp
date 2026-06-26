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

        Round_Down_Fractions(Xc, Yc, problem) ;
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
        xcost = Move_To_Vertex(x_flt, y_flt, gx, gy, Ai,Ap,c_flt,u,gamma,n,&j1);
        ycost = Move_To_Vertex(y_flt, x_flt, gy, gx, Ai,Ap,c_flt,u,gamma,n,&j2);

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


    /* ====================================================================== */
    // Phase 2: Greedy flips
    /* ====================================================================== */
    Round_Down_Fractions(X, Y, problem) ;
    Greedy_Flips(X, Y, problem, 9, (problem->level <= 2), false);


    /* ====================================================================== */
    // Phase 3: FM subspace gradient projection
    /* ====================================================================== */

#if 1

    /* Step 1: Convert x, y, c to floats, and compute gx and gy */
    x_flt = Point_To_Array(X->point) ;
    gx    = Az_To_Grad(c, gamma, Y->Az) ;

    y_flt = Point_To_Array(Y->point) ;    
    gy    = Az_To_Grad(c, gamma, X->Az) ;

    c_flt = Array_Int_To_Flt(c, n) ;

    /* Step 2: Do FM subspace method */
    // FM_Subspace()

    /* Step 3: Move x_flt, y_flt to almost binary points */
    int print_level = 0 ;

    VSPINT j1, j2 ;
    j1 = VSPboundary(x_flt, c_flt, gx, gy, n, Ai, Ap, gamma, print_level) ;
    j2 = VSPboundary(y_flt, c_flt, gy, gx, n, Ai, Ap, gamma, print_level) ;

    /* Step 4: Update the X and Y states with x_flt and y_flt */
    xcost = Move_To_Vertex(x_flt, y_flt, gx, gy, Ai, Ap, c_flt, u, gamma,n,&j1);
    ycost = Move_To_Vertex(y_flt, x_flt, gy, gx, Ai, Ap, c_flt, u, gamma,n,&j2);

    Array_To_State(x_flt, xcost, gy, X, problem, j1) ;
    Array_To_State(y_flt, ycost, gx, Y, problem, j2) ;

    /* Step 5: Cleanup */
    free(x_flt) ;
    free(y_flt) ;
    free(gx) ;
    free(gy) ;
    free(c_flt) ;
#endif

    /* ====================================================================== */
    // Phase 4: FM METIS
    /* ====================================================================== */
    n_passes = (l >= L / 2) ? 3 : 2 ;
    Round_Down_Fractions(X, Y, problem) ;

    FM_Metis(X, Y, problem, options, n_passes);
    if (Early_Exit(X, Y, problem)) goto CLEANUP;

    /* ====================================================================== */
    // Phase 5: Mountain climbing
    /* ====================================================================== */
    mountain_climbing(X, Y, problem, options, gamma, true) ;

    /* ====================================================================== */
    // Phase 6: Gamma reduction
    /* ====================================================================== */
    Gamma_Reduction_Loop(X, Y, problem, options);

    /* ====================================================================== */
    // Phase 7: Greedy flips
    /* ====================================================================== */
    Round_Down_Fractions(X, Y, problem) ;
    Greedy_Flips(X, Y, problem, 5, true, true);

    CLEANUP:

    #ifdef DEBUG
        check_states(X, Y, problem);
    #endif

    return;
}




// Only use this when the subproblem is much smaller than original problem
void VSP_FMsubspace
(
    VSPFLOAT  *x,    // variable in the bilinear objective
    VSPFLOAT  *y,    // variable in the bilinear objective
    VSPFLOAT  *a,    // vertex costs
    VSPINT    *Ai,   // array of graph neighbors
    VSPLONG   *Ap,   // column pointer for Ai
    VSPINT     u,    // upper bound for balance constraint
    VSPINT     n,    // number of vertices
    VSPFLOAT   gamma /*, max vertex cost
    VSPFLOAT   wx,      VSFLOAT workspace of size 7n
    VSPINT     wi,      VSINT workspace of size 3n + [nnz(A) +n]
    VSPLONG    wL       VSLINT workspace of size n+1 */
)
{
    int it, BigLoop, skip_napheapx, skip_napheapy, status ;
    VSINT Aip, i, j, k, improvement, nF, *F, *Si, *v, *wi ;
    VSLINT q, p, *Sp, *wL ;
    VSFLOAT eps, obj, obj_0, s, t, lambdax, lambday, ux, uy,
           *X, *Sx, *sx, *gx, *Gx, *Y, *Sy, *sy, *gy, *Gy, *gxtemp, *Sa, *wx ;

    /* The malloc's and their associated free's need to be moved to
       the main code. One malloc at the start of the code and one free
       at the completion. */
    wx = malloc (n*sizeof (VSPFLOAT)*12) ;
    wi = malloc (sizeof (VSINT)*Ap [n] + n*sizeof(VSINT)*3) ;
    wL = malloc ( ((VSLINT) n + 1)*sizeof(VSLINT)) ;

    ux = u ; // Initially, upper bound on x-balance constraint is u.
    uy = u ; /* Initially, upper bound on y-balance constraint is u.
                After the free variables are determined, the fixed variables
                (complement of the free variables) are removed from the
                problem and from u. Therefore ux and uy may be different,
                and both are generally smaller than u. */

    X  = wx     ; // X in subspace
    Sx = wx+1*n ; // free x variables forming subspace
    sx = wx+2*n ; // projection point in x
    gx = wx+3*n ; // total gradient wrt x in subspace
    Gx = wx+4*n ; // linear part of gradient wrt x in subspace
    Y  = wx+5*n ; // Y in subspace
    Sy = wx+6*n ; // free y variables forming subspace
    sy = wx+7*n ; // projection point in y
    gy = wx+8*n ; // total gradient wrt y in subspace
    Gy = wx+9*n ; // linear part of gradient wrt y in subspace
    gxtemp=wx+10*n;// gradient in x
    Sa = wx+11*n ; // a in subspace
    v = wi      ; // indices of free components
    F  = wi+1*n ; // F [i] >= 0 for free variables, F [i] = EMPTY for bound
    Si = wi+2*n ; // indices of free variables
    Sp = wL     ; // column pointers for Si

#if 1
VSFLOAT *GGX ;
GGX = malloc ( sizeof (VSFLOAT)*n ) ;
VSP_Grad (GGX, a, y, n, Ai, Ap, gamma, 1) ;
obj = VSP_Obj (x, GGX, y, a, n) ;
printf ("presubspace objective %e\n", obj) ;
free (GGX) ;
#endif

    BigLoop = 0 ;
    improvement = 1 ;
    vs_initi (F, EMPTY, n) ;
    //while ( improvement >= 1 )
    {
        BigLoop++ ;
        NAPdata *napdata ;
        napdata = napheap_setup () ;
    
        /* create subspace problem: vertices in separator and vertices directly
           connected to separator, includes fractional vertices */
        for (i = 0; i < n; i++)
        {
            // find vertices in separator
            if ( (x [i] == 0) && (y [i] == 0) )
            {
                F [i] = 0 ; // initialized as empty (diagonal is free)
                VSLINT const l = Ap [i+1] ;
                for (p = Ap [i]; p < l; p++)
                {
                    j = Ai [p] ;
                    // add to free set if not already free
                    if ( F [j] == EMPTY ) F [j] = 0 ;
                }
            }
        }
// !!! for testing, make all variables free !!!
//printf ("before printing matrix\n") ;
//fflush(stdout) ;
//for (i = 0; i < n; i++) F [i] = 0 ;
//vs_printAMATLAB (n, Ap, Ai, NULL, "A") ;
// printf ("after printing matrix\n") ;
//fflush(stdout) ;

        // nF counts number of free vertices
        nF = 0 ;
        // map from original vertices to subspace vertices
        for (i = 0; i < n; i++)
        {
            if ( F [i] != EMPTY ) // variable is free
            {
                F [i] = nF ; // maps original index to subspace
                v [nF] = i ; // maps subspace to original index
                Sa [nF] = a [i] ;
                Sx [nF] = x [i] ;
                Sy [nF] = y [i] ;
                nF++ ;
            }
            else // variable is fixed
            {
                ux -= x [i]*a [i] ;
                uy -= y [i]*a [i] ;
            }
        }
//printf ("nF: %i n: %i\n", nF, n) ;
//vs_printiMATLAB (F, n, "F") ;
//fflush(stdout) ;
     
        // build the subspace matrix and the gradients
        vs_initx (gx, VSZERO, nF) ;
        vs_initx (gy, VSZERO, nF) ;
        p = q = 0 ;
        Sp [p] = 0 ;
printf ("gamma: %e Sx [2] = %e\n", gamma, Sx [2]) ;
        for (j = 0; j < nF;)
        {
            s = t = VSZERO ; // obj deriv wrt yj or xj
            VSFLOAT const Sxj = -gamma*Sx [j] ; // xj in subspace times -gamma
            VSFLOAT const Syj = -gamma*Sy [j] ; // yj in subspace times -gamma
//printf ("Sxj: %e Syj: %e\n", Sxj, Syj) ;
//fflush(stdout) ;
            i = v [j] ;      // column of x or y in subspace matrix S
            VSLINT const l = Ap [i+1] ;
    
            if ( (p = Ap [i]) < l )
            {
                for (; p < l; p++)
                {
                    Aip = Ai [p] ;
//printf ("Aip: %i i: %i\n", Aip, i) ;
//fflush(stdout) ;
                    if ( Aip > i ) // add diagonal since not included in A
                    {
//printf ("add diagonal\n") ;
//fflush(stdout) ;
                        VSINT const k1 = F [i] ;
                        Si [q++] = k1 ;
                        if ( Syj != VSZERO ) gx [k1] += Syj ;
                        if ( Sxj != VSZERO ) gy [k1] += Sxj ;
                        VSINT const k2 = F [Aip] ;
                        if ( k2 != EMPTY )
                        {
                            Si [q++] = k2 ;
                            if ( Syj != VSZERO ) gx [k2] += Syj ;
                            if ( Sxj != VSZERO ) gy [k2] += Sxj ;
                        }
                        else // not free
                        {
                            t += y [Aip] ; // fixed part of x-gradient
                            s += x [Aip] ; // fixed part of y-gradient
if ( j == 2 ) printf ("add diagonal: Sxj = %e\n", gamma*x [Aip]) ;
                        }
                        p++ ;
                        break ;
                    }
                    VSINT const k = F [Aip] ;
                    if ( k != EMPTY )
                    {
                        Si [q++] = k ;
                        if ( Syj != VSZERO ) gx [k] += Syj ;
                        if ( Sxj != VSZERO ) gy [k] += Sxj ;
                    }
                    else // not free
                    {
                        t += y [Aip] ; // fixed part of x-gradient
                        s += x [Aip] ; // fixed part of y-gradient
if ( j == 2 ) printf ("not free1: Sxj = %e\n", gamma*x [Aip]) ;
                    }
                }
                for (; p < l; p++)
                {
                    Aip = Ai [p] ;
                    VSINT const k = F [Aip] ;
                    if ( k != EMPTY )
                    {
                        Si [q++] = k ;
                        if ( Syj != VSZERO ) gx [k] += Syj ;
                        if ( Sxj != VSZERO ) gy [k] += Sxj ;
                    }
                    else // not free
                    {
                        t += y [Aip] ; // fixed part of x-gradient
                        s += x [Aip] ; // fixed part of y-gradient
if ( j == 2 ) printf ("not free2: Sxj = %e\n", gamma*x [Aip]) ;
                    }
                }
                if ( Aip < i ) // add diagonal
                {
                    VSINT const k = F [i] ;
                    Si [q++] = k ;
                    if ( Syj != VSZERO ) gx [k] += Syj ;
                    if ( Sxj != VSZERO ) gy [k] += Sxj ;
                }
            }
            else // add diagonal
            {
                VSINT const k = F [i] ;
                Si [q++] = k ;
                if ( Syj != VSZERO ) gx [k] += Syj ;
                if ( Sxj != VSZERO ) gy [k] += Sxj ;
            }

            t = Sa [j] - gamma*t ;
            Gx [j] = t ;
            gx [j] += t ;
            s = Sa [j] - gamma*s ;
if ( j == 2 ) printf ("Sa [2] = %e\n", Sa [j]) ;
            Gy [j] = s ;
            gy [j] += s ;
            j++ ;
            Sp [j] = q ;
        }
#if 0
if ( Sp [n] != Ap [n]+n )
{
    printf ("error Sp [%i] = %i != Ap [%i] = %i\n", n, Sp [n], n, Ap [n]) ;
    fflush(stdout) ;
    exit(1) ;
}
q = p = 0 ;
for (i = 0; i < n; i++)
{
    VSLINT const l = Sp [i+1] ;
    for (; p < l; p++)
    {
        if ( Si [p] == i ) continue ;
        if ( Si [p] != Ai [q] )
        {
            printf ("error at row: %i Si [%i] = %i != Ai [%i] = %i\n",
                     i, p, Si [p], q, Ai [q]) ;
            exit(1) ;
        }
        q++ ;
    }
}
for (i = 0; i <= n; i++)
{
    if ( Ap [i]+i != Sp [i] )
    {
        printf ("error at row: %i Ap [%i] + %i = %i != Sp [%i] = %i\n",
                 i, i, i, Ap [i], i, Sp [i]) ;
        exit(1) ;
    }
}
        napdata->a = NULL ;
        for (i = 0; i < n; i++)
        {
            if ( a [i] != 1 )
            {
                napdata->a = a ;
                break ;
            }
        }
#endif
#if 1
// check that Si and Sp includes the diagonal; indices in column increasing
p = 0 ;
for (j = 0; j < nF; j++)
{
    VSLINT const l = Sp [j+1] ;
    i = -1 ;
    k = 0 ;
    for (; p < l; p++)
    {
        if ( Si [p] <= i ) // i stores last index in column
        {
            printf ("indices do not increase at col = %i p = %i\n", j, p) ;
            fflush (stdout) ;
            exit (1) ;
        }
        i = Si [p] ;
        if ( i == j ) k = 1 ; // diagonal found
    }
    if ( k == 0 )
    {
        printf ("diagonal not found in col = %i\n", j) ;
        fflush (stdout) ;
        exit (1) ;
    }
}

// check Gx, Gy, Si, and Sp
VSFLOAT *GGx, *GGy ;
GGx = malloc ( sizeof (VSFLOAT)*nF ) ;
vs_copyx (GGx, Sa, nF) ;
GGy = malloc ( sizeof (VSFLOAT)*nF ) ;
vs_copyx (GGy, Sa, nF) ;
p = 0 ;
for (j = 0; j < n; j++)
{
    VSLINT const l = Ap [j+1] ;
    int ifree, jfree, nfree ;
    nfree = 0 ;
    jfree = -1 ;
    if ( F [j] != -1 )
    {
        jfree = F [j] ;
        q = Sp [jfree] ;
    }
    for (; p < l; p++)
    {
        ifree = -1 ;
        i = Ai [p] ;
        if ( F [i] != -1 ) ifree = F [i] ;
        if ( (ifree != -1) && (jfree != -1) ) // in Si and Sp
        {
            if ( q >= Sp [jfree+1] )
            {
                printf ("reached end of col %i in S without finding %i\n",
                         jfree, ifree) ;
                printf ("which corresponds to (%i, %i) in A\n", i, j) ;
                fflush (stdout) ;
                exit(1) ;
            }
            nfree++ ;
            // find row ifree in col jfree
            if ( Si [q] == jfree ) q++ ; // skip diagonal of S
            if ( Si [q] != ifree )
            {
                    printf ("A(%i, %i) not in col %i of S\n", i, j, jfree) ;
                    fflush (stdout) ;
                    exit (1) ;
            }
            q++ ;
        }
        if ( (ifree != -1) && (jfree == -1) ) // belongs in GGx
        {
            GGx [ifree] -= gamma*y [j] ;
        }
        if ( (ifree == -1) && (jfree != -1) ) // belongs in GGy
        {
            GGy [jfree] -= gamma*x [i] ;
        }
    }
    if ( (jfree != EMPTY) && Sp [jfree+1] - Sp [jfree] != nfree + 1 )
    {
        printf ("col %i of A has %i free elements while S "
                "has %i elements\n",
                 j, nfree, Sp [jfree+1]-Sp [jfree] - 1) ;
        fflush (stdout) ;
        exit (1) ;
    }
}
// check GGx and GGy
s = t = 0 ;
for (i = 0; i < nF; i++)
{
    t += fabs (GGx [i] - Gx [i]) ;
    s += fabs (Gx [i]) ;
}
printf ("relative error %e in Gx\n", t/s) ;
if (t/s >= 1.e-10 )
{
    printf ("relative error %e in Gx\n", t/s) ;
    fflush (stdout) ;
    exit (1) ;
}
s = t = 0 ;
for (i = 0; i < nF; i++)
{
    t += fabs (GGy [i] - Gy [i]) ;
    s += fabs (Gy [i]) ;
}
printf ("relative error %e in Gy\n", t/s) ;
if (t/s >= 1.e-10 )
{
    printf ("relative error %e in Gy\n", t/s) ;
    fflush (stdout) ;
    exit (1) ;
}
free (GGx) ;
free (GGy) ;
#endif
        napdata->n = nF ;   // problem dimension at current level
        napdata->a = Sa ;   // linear constraint for subspace
        napdata->blo = 0 ;  // lower bound of balance constraint
        eps = nF*1.e-8 ;    // when cost does not improve by eps, exit
    
        /* initial objective, obj. This is the objective for the subproblem
           with the constant term in the original objective dropped */
        obj = VSP_Obj (Sx, gx, Sy, Gy, nF) ;
printf ("gamma = %e\n", gamma) ;
#if 0
printf ("nF: %i gamma: %e obj: %e\n", nF, gamma, obj) ;
vs_printxMATLAB(Sx, nF, "x") ;
vs_printxMATLAB(Sy, nF, "y") ;
vs_printxMATLAB(Sa, nF, "a") ;
vs_printxMATLAB(gx, nF, "gx") ;
for (i = 0; i < n; i++) printf ("F[%i] = %i\n", i, F [i]) ;
for (i = 0; i < n; i++) printf ("v[%i] = %i\n", i, v [i]) ;
fflush(stdout) ;
#endif
        VSFLOAT start_obj = obj ;
    
        // starting guess for lambda
        lambdax = VSZERO ;
        lambday = VSZERO ;
    
        /* Perform iterations of the gradient projection algorithm for subspace
           The new (Sx, Sy) is the solution of
           
                  min 0.5*(x-sx)^2 + 0.5*(y-sy)%2
          
                       subject to 0 <= ax'x <= ux, 0 <= x <= 1,
                                       ay'y <= uy, 0 <= y <= 1,
        
           where sx and sy are steps along the gradient starting from the
           iterates Sx and Sy: sx = Sx + s*gx and sy = Sy + s*gy,
           s = stepsize */

        // check that gx and gy are gradients in the subspace
GGx = malloc ( sizeof (VSFLOAT)*nF ) ;
GGy = malloc ( sizeof (VSFLOAT)*nF ) ;
VSP_Grad (GGx, Gx, Sy, nF, Si, Sp, gamma, 0) ;
s = t = 0 ;
for (i = 0; i < nF; i++)
{
    t += fabs (GGx [i] - gx [i]) ;
    s += fabs (gx [i]) ;
}
printf ("relative error %e in gx\n", t/s) ;
if (t/s >= 1.e-10 )
{
    printf ("relative error %e in gx\n", t/s) ;
    fflush (stdout) ;
    exit (1) ;
}
VSP_Grad (GGy, Gy, Sx, nF, Si, Sp, gamma, 0) ;
s = t = 0 ;
for (i = 0; i < nF; i++)
{
    t += fabs (GGy [i] - gy [i]) ;
    s += fabs (gy [i]) ;
}
printf ("relative error %e in gy\n", t/s) ;
if (t/s >= 1.e-10 )
{
    printf ("relative error %e in gy\n", t/s) ;
    fflush (stdout) ;
    exit (1) ;
}
free (GGx) ;
free (GGy) ;
    
        // start of the gradient projection algorithm loop
        it = 0 ;
        while ( 1 )
        {
            obj_0 = obj ; // save prior objective value
            s = vs_max(gx, nF) ;
            t = vs_max(gy, nF) ;
            s = MAX(s, t) ;
            // if gradient vanishes, then current point is stationary
            if ( s == VSZERO ) break ;
            it++ ;
    
            // Scaling factor for gradients
            s = 0.8909725/s ;
            // sfac = factor by which stepsize decreases to improve obj
            VSFLOAT sfac = 0.2 ;
    
            /* Project gradient onto feasible set. The projection decomposes
               into two separate projections: sx is projected on the x-feasible
               set and sy is projected on the y-feasible set.  */
            obj = -NAPINF ; /* initialize obj to a lower bound on objective*/
            int subit = 0 ;
#if 0
vs_printAMATLAB (nF, Sp, Si, NULL, "S") ;
vs_printxMATLAB (Gx, nF, "Gx") ;
vs_printxMATLAB (Gy, nF, "Gy") ;
vs_printxMATLAB (gx, nF, "gx") ;
vs_printxMATLAB (gy, nF, "gy") ;
vs_printxMATLAB (Sx, nF, "x") ;
vs_printxMATLAB (Sy, nF, "y") ;
#endif
            while(obj + eps < obj_0 )
            {
                subit++ ; 
                if ( subit == 1 )
                {
                    // scale gradient by optimal stepsize 0.891/length
                    vs_scale (gx, s, nF) ;
                    vs_scale (gy, s, nF) ;
                }
                else
                {
                    vs_scale (gx, sfac, nF) ;
                    vs_scale (gy, sfac, nF) ;
                }
    
                // Compute new point (sx, sy)
                vs_vecadd(sx, Sx, gx, nF) ;
                vs_vecadd(sy, Sy, gy, nF) ;
#if 0
vs_printxMATLAB (Sx, nF, "x1") ;
vs_printxMATLAB (gx, nF, "gx1") ;
vs_printxMATLAB (sx, nF, "sx1") ;
vs_printxMATLAB (Sy, n, "y1") ;
vs_printxMATLAB (gy, n, "gy1") ;
vs_printxMATLAB (sy, n, "sy1") ;
#endif
    
#if 0
printf ("new point:\n") ;
vs_printxMATLAB (sx, nF, "sx") ;
vs_printxMATLAB (sy, nF, "sy") ;
printf ("lambdax: %e\n", lambdax) ;
#endif
                // Project sx to feasible set, see if lambda = 0 works
                skip_napheapx = VSFALSE ;
                if ( !lambdax )
                {
                    t = VSZERO ;
                    for (i = 0; i < nF; i++)
                    {
                        if ( sx [i] > VSONE )
                        {
                            sx [i] = VSONE ;
                            t += Sa [i] ;
                        }
                        else if ( sx [i] < VSZERO )
                        {
                            sx [i] = VSZERO ;
                        }
                        else
                        {
                            t += Sa [i]*sx [i] ;
                        }
                    }
                    if ( t <= ux )
                    {
                        skip_napheapx = VSTRUE ;
                    }
                    else
                    {
                        // restore sx to correct value
                        vs_vecadd(sx, Sx, gx, nF) ;
                    }
                }
                if ( !skip_napheapx )
                {
//printf ("use napheap x\n") ;
                    napdata->c = sx ;   // linear term in projection
                    napdata->x = X ;    // solution returned in X
                    napdata->bhi = ux ; // upper bound of balance constraint
    
                    napdata->lambda = lambdax ;
                    status = napheap(napdata) ;
                    lambdax = napdata->lambda ;
#if 0
                    //printf ("lambdax: %e\n", lambdax) ;
                    //if( status ) napheap_print_status (napdata) ;
                    VSFLOAT err, r ;            
                    err = napheap_check (&r, &s, &t, napdata) ;
                    printf ("errx: %e errb: %e errB: %e errc: %e\n",
                             err, r, s, t) ;
                    if ( err > 1.e-6 )
                    {
                        printf ("NAPHEAP error: %e\n", err) ;
                        exit(1) ;
                    }
#endif
                }
    
                // Project sy to feasible set, see if lambda = 0 works
                skip_napheapy = VSFALSE ;
//printf ("lambday: %e\n", lambday) ;
                if ( !lambday )
                {
                    t = VSZERO ;
                    for (i = 0; i < nF; i++)
                    {
                        if ( sy [i] > VSONE )
                        {
                            sy [i] = VSONE ;
                            t += Sa [i] ;
                        }
                        else if ( sy [i] < VSZERO )
                        {
                            sy [i] = VSZERO ;
                        }
                        else
                        {
                            t += Sa [i]*sy [i] ;
                        }
                    }
                    if ( t <= uy )
                    {
                        skip_napheapy = VSTRUE ;
                    }
                    else
                    {
                        // restore sy to correct value
                        vs_vecadd(sy, Sy, gy, nF) ;
                    }
                }
                if ( !skip_napheapy )
                {
//printf ("use napheap y\n") ;
                    napdata->c = sy ;   // linear term in projection
                    napdata->x = Y ;    // solution returned in Y
                    napdata->bhi = uy ; // upper bound of balance constraint
                    napdata->lambda = lambday ;
                    status = napheap(napdata) ;
                    lambday = napdata->lambda ;
#if 0
                    //printf ("lambday: %e\n", lambday) ;
                    //printf ("n = %i ;\n", n) ;
                    //printf ("u = %i ;\n", u) ;
                    //printf ("l = 0 ;\n") ;
                    //printf ("lambda after = %25.15e\n", napdata->lambda) ;
                    //int i ;
                    //for (i = 0; i < n; i++)
                    //{
                        //    printf ("y[%i] = %25.15e ;\n", i, sy [i]) ;
                    //}
                    //for (i = 0; i < n; i++)
                    //{
                    //    printf ("a[%i] = %25.15e ;\n", i, a [i]) ;
                    //}
                    //if( status ) napheap_print_status (napdata) ;
    
                    err = napheap_check (&r, &s, &t, napdata) ;
                    printf ("erry: %e errb: %e errB: %e errc: %e\n",
                             err, r, s, t) ;
                    if ( err > 1.e-6 )
                    {
                        printf ("NAPHEAP error: %e\n", err) ;
                        exit(1) ;
                    }
#endif
                }
    
                // compute gxtemp
                if ( skip_napheapy )
                {
                    // solution stored in sy
                    VSP_Grad (gxtemp, Gx, sy, nF, Si, Sp, gamma, 0) ;
                    if ( skip_napheapx )
                    {
                        // solution stored in sx
                        obj = VSP_Obj (sx, gxtemp, sy, Gy, nF) ;
                    }
                    else
                    {
                        // solution stored in X
                        obj = VSP_Obj (X, gxtemp, sy, Gy, nF) ;
                    }
                }
                else
                {
                    // solution stored in Y
                    VSP_Grad (gxtemp, Gx, Y, nF, Si, Sp, gamma, 0) ;
                    if ( skip_napheapx )
                    {
                        // solution stored in sx
                        obj = VSP_Obj (sx, gxtemp, Y, Gy, nF) ;
                    }
                    else
                    {
                        // solution stored in X
                        obj = VSP_Obj (X, gxtemp, Y, Gy, nF) ;
                    }
                }
    
                //printf("it: %i obj: %f    obj0: %f\n", it, obj, obj_0) ;
                printf ("it: %i obj0: %f obj: %f +: %f n: %i nF: %i\n",
                           it, obj_0, obj, obj-obj_0, n, nF);
//if ( obj < obj_0 )
//{
    //fflush(stdout) ;
    //exit(1) ;
//}
            }
    
            /* Update: x <- X and y <- Y */
            if ( skip_napheapx )
            {
                // solution stored in sx
                vs_copyx(Sx, sx, nF) ;
            }
            else
            {
                // solution stored in X
                vs_copyx(Sx, X, nF) ;
            }
            if ( skip_napheapy )
            {
                // solution stored in sy
                vs_copyx(Sy, sy, nF) ;
            }
            else
            {
                // solution stored in Y
                vs_copyx(Sy, Y, nF) ;
            }
    
            vs_copyx(gx, gxtemp, nF) ;        // grad at x stored in gxtemp
    
            // evaluate gradient at y
            VSP_Grad (gy, Gy, Sx, nF, Si, Sp, gamma, 0) ;
            // |obj - obj_0| < eps => stationary
            if ( obj < obj_0 + eps ) break ;
            if ( obj - obj_0 <= 0.005*VSMAX(obj, 1) ) break ;
        }
        printf ("start: %f  end: %f +: %f\n",
                 start_obj, obj, obj - start_obj) ;
        improvement = obj - start_obj ;
#if 0
        for (j = 0; j < nF; j++)
        {
            i = v [j] ; // column for x or y
            x [i] = Sx [i] ;
            y [i] = Sy [i] ;
            F [i] = EMPTY ;
        }
        
GGx = malloc ( sizeof (VSFLOAT)*n ) ;
VSP_Grad (GGx, a, y, n, Ai, Ap, gamma, 1) ;
obj = VSP_Obj (x, GGx, y, a, n) ;
free (GGx) ;
printf ("postsubspace objective %e\n", obj) ;
fflush(stdout) ;
exit(1) ;
#endif

        // convert Sx to almost binary
        VSPboundary (Sx, Sa, gx, nF, 0) ;
        // convert Sy to almost binary after updating y's gradient
        VSP_Grad (gy, Gy, Sx, nF, Si, Sp, gamma, 0) ;
        VSPboundary (Sy, Sa, gy, nF, 0) ;
 
        // copy Sx and Sy back into original solution x
        for (j = 0; j < nF; j++)
        {
            i = v [j] ; // column for x or y
            x [i] = Sx [i] ;
            y [i] = Sy [i] ;
            F [i] = EMPTY ;
        }
        
    
        napheap_terminate(&napdata) ;
    }
//printf ("end of routine\n") ;
//fflush(stdout) ;
    free (wx) ;
    free (wi) ;
    free (wL) ;
}

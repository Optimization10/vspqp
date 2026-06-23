/* =========================================================================
   ============================= vs (VSP_STARTUP) ==========================
   =========================================================================
       ________________________________________________________________
      |    Generate a starting guess for the vertex separator problem  |
      |       from an approximation to the solution of the graph       |
      |       bisection problem. The graph bisection problem can be    |
      |       formulated as                                            |
      |                                                                |
      |                   1                                            |
      |    (GB)       min - (1-x)'Ax subject to a'x = m, 0 <= x <= 1   |
      |                   2                                            |
      |                                                                |
      |    n = number of vertices in an undirected graph,              |
      |    m = n/2, elements of "a" vector are positive and sum to n.  |
      |    A = D + E, E is the edge weighted adjacency matrix for graph|
      |    D = diagonal matrix with d_{jj} maximum element in column j |
      |    of E. The constraint 0 <= x <= 1 is replaced by ||x-c|| <= r|
      |    where c is the center of the box, the point whose components|
      |    are all 1/2. After the change of variables x = y + c, this  | 
      |    sphere constrained problem is                               |
      |                                                                |
      |                     1                                          |
      |    (y)        min - -  y'Ay subject to a'y = 0, ||y|| <= r.    |
      |                     2                                          |
      |                                                                |
      |    This problem is equivalent to                               |
      |                                                                |
      |                     1                                          |
      |    (z)        min - -  (Lz)'A(Lz) subject to ||Lz|| <= r       |
      |                     2                                          |
      |                                                                |
      |    where L is the matrix that projects a vector into the space |
      |    orthogonal to "a". The vector z = a is an eigenvector of LAL|
      |    corresponding to the eigenvalue 0. A solution of (y) is any |
      |    eigenvector LAL orthogonal to "a" associated with the       |
      |    smallest eigenvalue of LAL. We approximate this eigenvector |
      |    using Lanczos method. The corresponding approximation to a  |
      |    solution of (GB) is x_i = 0 (z_i < 0) and x_i = 1 (z_i > 0) |
      |    If z_i = 0, then vertex i is put in the separator. We treat |
      |    S = { i : x_i = 0 } and T = { i : x_i = 1 } as a partition  |
      |    of the vertices of the graph. The vertex set S_c subset of S|
      |    and T_c of subset T are the set of vertices associated with |
      |    cut edges. That is, for each i in S_c, there exist j in T_c |
      |    such that e_{ij} != 0. The Ashcraft/Liu algorithm is used to|
      |    find a minimum vertex cover. These vertices form the        |
      |    separator set together with the vertices where z_i = 0      |
      |                                                                |
      |                     Version 1.0 (January 1, 2025)              |
      |                                                                |
      |                         William W. Hager                       |
      |                           hager@ufl.edu                        |
      |                          Ovidiu Nechita                        |
      |                         onechita@ufl.edu                       |
      |                     Department of Mathematics                  |
      |                       University of Florida                    |
      |                     Gainesville, Florida 32611                 |
      |                                                                |
      |                   Copyright by William W. Hager                |
      |________________________________________________________________|

       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|*/

#include "vs_internal.h"
#include <math.h>
#include <stdio.h>

int vs
(
    int       *x, /* x [i] = 0, 1, or 2 to put i in set 0, set 1, or separator
                     respectively,  malloc'd by user. For VSP-based start
                     guess, x is starting guess for bisection problem.
                     If x [i] = 2, then use 1/2 for the associated starting
                     guess for bisection problem. */
    VSINT    *Ei, /* row numbers */
    VSLINT   *Ep, /* column pointers */
    VSINT    *Ex, /* edge weights,   NULL => all entries in Ex are 1 */
    VSINT    *vw, /* vertex weights, NULL => all entries in a are 1 */
    VSINT      n, /* problem dimension */
    VSFLOAT    m, /* balance constraint parameter: 0.5 * original # vertices
                     m = 0.5*sum(vw [i]) if vw is provided */
    VSparm *Parm  /* parameter structure */
)
{
    int qr_its ;
    VScom Com ;
    Com.Parm = Parm ;

    VSINT n0, n1, ns, cutedge ;
    VSINT i, j, j0, *Ai ;
    VSLINT nz, nzA, p, q, *Ap ;
    VSFLOAT Amax, r, s, t, *a, *A12, *Ax, *b, *z ;

    int status = 0 ;
    int const PrintLevel = Parm->PrintLevel ;
    if ( (PrintLevel > 0) && Parm->PrintParms ) vs_print_parms (Parm) ;

    int const vw_is_one = (vw  == NULL) ? VSTRUE : VSFALSE ;
    int const E_is_one = (Ex == NULL) ? VSTRUE : VSFALSE ;

    /* Ax is only required if USE_MKL is true or if Ex is not NULL */
    int need0Ax = !E_is_one ;
#ifdef USE_MKL
    need0Ax = VSTRUE ;
#endif
    int const needAx = need0Ax ;

    /* convert vertex weight array to VSFLOAT array if vw not NULL */
    if ( vw_is_one ) a = NULL ;
    else             a = vs_malloc (&status, n, sizeof (VSFLOAT)) ;

    /* incorporate diagonal in A, nz = number of edges in graph */
    nz  = Ep [n] ;
    /* nz's in A including diagonal */
    nzA = nz + n ;
    /* row indices */
    Ai = vs_lmalloc (&status, nzA, sizeof (VSINT)) ;
    Ax = NULL ;
    if ( needAx )
    {
        Ax = vs_lmalloc (&status, nzA, sizeof (VSFLOAT)) ;
    }
    /* row/column pointers */
    Ap = vs_malloc (&status, n+1, sizeof (VSLINT)) ;/* pointers */
    /* z = solution of Lz prob */
    z  = vs_malloc (&status, n, sizeof(VSFLOAT)) ;
    /* return if out-of-memory */
    if ( status ) return (status) ;
    if ( !vw_is_one ) for (i = 0; i < n; i++) a [i] = vw [i] ;
    VSFLOAT const aa = (vw_is_one) ? VSONE/n : VSONE/vs_dot (a, a, n) ;
    p = 0 ; /* points into E */
    q = 0 ; /* points into A */
    Amax = 0 ; /* largest element in A */
    j0 = 0 ; /* first column to process */
    Ap [0] = 0 ;
    /* build Ax which is the same as Ex except for diagonal element =
       largest element of E in column */
    for (j = 1; j <= n; )
    {
        VSLINT const l = Ep [j] ;
        if ( !needAx ) /* do not need Ax */
        {
            Amax = 1 ;
            VSLINT const qstart = q ;
            for (; p < l; )
            {
                if ( Ei [p] > j0 )    /* beyond diagonal of E */
                {
                    Ai [q++] = j0 ;
                    break ;
                }
                else
                {
                    Ai [q++] = Ei [p++] ;
                }
            }
            if ( p == l ) /* no break occurred, insert diagonal element */
            {
                Ai [q++] = j0 ;
            }
            else          /* break occurred and diagonal was added */
            {
                for (; p < l; )
                {
                    Ai [q++] = Ei [p++] ;
                }
            }
            Ap [j] = q ;
            j0 = j++ ;
        }
        else /* need Ax */
        {
            /* stores largest element in column */
            if ( E_is_one && (p < l) ) t = VSONE ;
            else                       t = VSZERO ;
            
            for (; p < l; )
            {
                /* beyond diagonal of E, will break loop */
                if ( Ei [p] > j0 )
                {
                    i = q ; /* fill in Ax [i] at end of loop */
                    Ai [q++] = j0 ;
                    break ;
                }
                else
                {
                    if ( E_is_one )
                    {
                        Ax [q] = -VSONE ;
                    }
                    else
                    {
                        if ( Ex [p] > t ) t = Ex [p] ;
                        Ax [q] = -Ex [p] ;
                    }
                    Ai [q++] = Ei [p++] ;
                }
            }
            if ( p == l ) /* no break occurred, insert diagonal */
            {
                Ax [q] = -t ;
                if ( t > Amax ) Amax = t ;
                Ai [q++] = j0 ;
            }
            else
            {
                for (; p < l; )
                {
                    if ( E_is_one ) Ax [q] = -VSONE ;
                    else
                    {
                        if ( Ex [p] > t ) t = Ex [p] ;
                        Ax [q] = -Ex [p] ;
                    }
                    Ai [q++] = Ei [p++] ;
                }
                Ax [i] = -t ;
                if ( t > Amax ) Amax = t ;
            }
            Ap [j] = q ;
            j0 = j++ ;
        }
    }

    /* Two different approaches are used to compute a starting guess for
       the VSP. Approach 1, as explained in the documentation,
       involves solving an eigenvalue problem associated with the bisection
       problem. This approach is used at a coarse graph, before we have an
       approximation to the solution of the VSP. In Approach 2,
       we use the starting guess for the VSP, which was computed by
       interpolation from a coarser graph, to generate a starting guess
       for the bisection problem. Only use this approach when one of the 
       variables associated with the VSP, either x or y, is near an upper
       or lower bound of the balance constraint. When a variable is near
       such a bound, this bound reduces the search space.
       This unbalanced starting guess for the bisection problem leads to
       a balanced solution of the bisection problem, which yields
       a balanced starting point for the VSP, which is less restricted by
       the balance constraint.

       The next block of code is used to compute the eigenvector-based
       starting guess for the bisection problem. */
    b = (VSFLOAT *) vs_malloc (&status, n, sizeof (VSFLOAT)) ;
    Com.wx = b ;
    Com.wi = (VSINT *)   vs_malloc (&status, n, sizeof (VSINT)) ;
    if ( Parm->compute_eig )
    {
        VSINT LO, nO, nOmax, nQR, nRot, SmallMatrix ;
        VSFLOAT ab, ba, Lanczos_tol,
               *A, *AP, *d, *u, *V, *vnew, *vprior, *vpprior ;
    

        VSgivens Gstruc ;
        VSgivens *G = &Gstruc ;
        Com.G = G ;

        /* qr_its times matrix dimension = max # QR iter */
        qr_its = Parm->qr_its ;
        SmallMatrix = Parm->SmallMatrix ;
        /* allocate memory for work arrays */
        if ( n <= SmallMatrix ) /* small matrix */
        {
            /* # columns in orthogonal basis, nO = actual # cols,
                                           nOmax = max cols*/
            nO = nOmax = n ;
        }
        else /* large matrix, use Lanczos process */
        {
            int n2 ;
            /* n2 = log2(n) */
            frexp((double) n, &n2) ;
            /* max # orthogonal basis vectors */
            nOmax = Parm->basis_factor*(n2-1) ;
            /* actual number nO of basis vector <= nOmax */
            V = (VSFLOAT *) vs_Lmalloc (&status, (nOmax+1)*(int64_t) n,
                                         sizeof (VSFLOAT)) ;
        }
        /* max number of QR iterations */
        nQR = nOmax*qr_its ;   /* QR iter = max dimension * qr_its */
        nRot = nOmax*nQR ;     /* max # Givens rotations */
    
        /* space for final dense matrix */
        A = (VSFLOAT *) vs_malloc(&status, nOmax*nOmax, sizeof (VSFLOAT));
    
        /* allocate space for Givens rotations in diagonalization */
        G->gx = vs_malloc (&status, nRot, sizeof (VSFLOAT)) ;
        G->gy = vs_malloc (&status, nRot, sizeof (VSFLOAT)) ;
        G->gz = vs_malloc (&status, nRot, sizeof (short)) ;
        G->gs = vs_malloc (&status, (nRot/30)+nQR+nOmax+5, sizeof (VSFLOAT)) ;
        G->gi = vs_malloc (&status, (nRot/30)+nQR+nOmax+5, sizeof (VSINT)) ;
        G->gj = vs_malloc (&status, nQR, sizeof (VSINT)) ;
        G->gk = vs_malloc (&status, nQR, sizeof (VSINT)) ;
        G->gf = vs_malloc (&status, nQR, sizeof (VSINT)) ;
        d = (VSFLOAT *) vs_malloc (&status, nOmax, sizeof (VSFLOAT)) ;
        u = (VSFLOAT *) vs_malloc (&status, nOmax, sizeof (VSFLOAT)) ;
        if ( status ) return (status) ;
    
        /* tridiagonalize the matrix */
        if ( n <= SmallMatrix )
        {
            /* Small matrices are tridiagonalized by Householder reflections.
               First, create the dense matrix
    
                   LAL = [I - aa'/(a'a)] * A * [I - aa'/(a'a)]
    
               If b = A*a, then
    
                   LAL = A + a[a(b'a)/(a'a)^2 - b/(a'a)]' - ba'/(a'a) ,
    
               which can be computed column by column. */
    
            /* evaluate b and ba := b'a */
    
            ba = VSZERO ;
            if ( E_is_one && vw_is_one )
            {
                for (j = 0; j < n; j++)
                {
                    t = Ap [j] - Ap [j+1] ; /* < 0 */
                    b [j] = t ;
                    ba += t ;
                }
            }
            else if ( vw_is_one ) /* vw = 1, but E not 1 */
            {
                p = 0 ;
                for (j = 0; j < n; j++)
                {
                    VSLINT const l = Ap [j+1] ;
                    t = VSZERO ;
                    for (; p < l; p++)
                    {
                        t += Ax [p] ;
                    }
                    b [j] = t ;
                    ba += t ;
                }
            }
            else if ( !vw_is_one && E_is_one ) /* vw_not_one && E_is_one */
            {
                p = 0 ;
                for (j = 0; j < n; j++)
                {
                    VSLINT const l = Ap [j+1] ;
                    t = VSZERO ;
                    for (; p < l; p++)
                    {
                        t -= a [Ai [p]] ;
                    }
                    b [j] = t ;
                    ba += t*a [j] ;
                }
            }
            else /* general case: both vw and E not one */
            {
                p = 0 ;
                for (j = 0; j < n; j++)
                {
                    VSLINT const l = Ap [j+1] ;
                    t = VSZERO ;
                    for (; p < l; p++)
                    {
                        t += Ax [p]*a [Ai [p]] ;
                    }
                    b [j] = t ;
                    ba += t*a [j] ;
                }
            }
    
            /* construct the dense A matrix column by columns */
            t = ba*aa*aa ; /* b'a/(a'a)^2 since aa = 1/(a'a) */
            AP = A ;
            p = 0 ;
            /* LAL = A + a[a(b'a)/(a'a)^2 - b/(a'a)]' - ba'/(a'a) */
            for (j = 0; j < n;)
            {
                if ( vw_is_one )
                {
                    s = t - b [j]*aa ;
                    for (i = 0; i < n; i++) AP [i] = s - b [i]*aa ;
                    }
                else /* !vw_is_one */
                {
                    s = a [j]*t - b [j]*aa ;
                    r = a [j]*aa ;
                    for (i = 0; i < n; i++) AP [i] = a [i]*s - b [i]*r ;
                }
                /* subtract Ax from current column */
                VSLINT const l = Ap [++j] ;
    
                if ( E_is_one )
                {
                    for (; p < l; p++) AP [Ai [p]]-- ;
                }
                else
                {
                    for (; p < l; p++) AP [Ai [p]] += Ax [p] ;
                }
    
                AP += n ;
            }
    
            /* reduce A to tridiagonal form using Householder matrices */
            vs_house (A, d, u, n, &Com) ;
    
            /* diagonalize the tridiagonal matrix using Givens rotations */
            vs_diag (d, u, n, &Com) ;
    
            /* find the smallest eigenvalue */
            t = d [0] ;
            i = 0 ;
            for (j = 1; j < n; j++)
            {
                if ( d [j] < t )
                {
                    t = d [j] ;
                    i = j ;
                }
            }
            /* t = smallest eigenvalue of L'AL and i = index in d,
               evaluate corresponding eigenvector of L'AL */
            vs_initx (b, VSZERO, n) ;
            b [i] = VSONE ;
            vs_tGivensMult (b, Com.G, Com.npass, Com.nscale, n) ;
            vs_DenseMult (z, b, A, n, n) ;
        }
        else /* n > SmallMatrix, Lanczos tridiagonalization */
        {
            /* stop Lanczos iteration when |u [j]| <= Amax*(machine epsilon) */
            Lanczos_tol = Amax*Parm->Lanczos_tol ;
    
            /* seed for random number generator, default is 1 */
            srand (Parm->seed) ;
            VSFLOAT rand_norm = VSONE/RAND_MAX ;
    
            /* random starting guess:
               first choose random numbers on [-0.5, +0.5] */
            for (j = 0; j < n; j++)
            {
                V [j] = rand_norm*rand() - 0.5 ;
            }
            /* replace starting vector by it projection */
            vs_project (V, a, aa, vw_is_one, n) ;
            t = vs_norm (V, n) ;
            t = VSONE/t ;
            /* normalize the starting vector */
            vs_scale (V, t, n) ;
            vnew = V+n ;
            /* vprior = random starting guess projected on orth complement (V)*/
            vprior = V ;
            /* Lanczos startup process: vnew = LA*vprior */
            vs_Amult (vnew, vprior, Ap, Ai, Ax, n) ;
            vs_project (vnew, a, aa, vw_is_one, n) ; /* vnew = LAL * vprior */
            d [0] = vs_dot (vnew, vprior, n) ;
            vs_daxpy (vnew, vprior, -d [0], n) ;
            u [0] = vs_norm (vnew, n) ;
            /* evaluate second Lanczos vector if termination not reached */
            if ( fabs(u [0]) > Lanczos_tol )
            {
                /* normalize vnew to be a unit vector */
                t = VSONE/u [0] ;
                vs_scale (vnew, t, n) ;
                nO = 1 ;
                while ( 1 )
                {
                    vpprior = vprior ;
                    vprior = vnew ;
                    vnew = vnew+n ;
                    /* vnew = LAL * vprior where L*vprior = vprior
                       initially compute vnew = A*vprior */
                    vs_Amult (vnew, vprior, Ap, Ai, Ax, n) ;
                    /* now project vnew to obtain vnew = L*A*vprior */
                    vs_project (vnew, a, aa, vw_is_one, n) ;
                    d [nO] = vs_dot (vnew, vprior, n) ;
                    nO++ ;
                    if ( nO >= nOmax ) break ;
                    vs_daxpy (vnew, vpprior, -u [nO-2], n) ;
                    vs_daxpy (vnew, vprior, -d [nO-1], n) ;
                    s = u [nO-1] = vs_norm (vnew, n) ;
                    if ( s <= Lanczos_tol ) break ;
                    s = VSONE/s ;
                    vs_scale (vnew, s, n) ;
                }
            }

            /* diagonalize the tridiagonal matrix using Givens rotations */
            vs_diag (d, u, nO, &Com) ;
    
            /* find the smallest eigenvalue */
            t = d [0] ;
            i = 0 ;
            for (j = 1; j < nO; j++)
            {
                if ( d [j] < t )
                {
                    t = d [j] ;
                    i = j ;
                }
            }
            /*printf ("eval = %25.15e index = %i\n", t, i) ;*/
            /* t = smallest eigenvalue of L'AL and i = index in d,
               evaluate corresponding eigenvector of L'AL */
            vs_initx (b, VSZERO, nO) ;
            b [i] = VSONE ;
            vs_tGivensMult (b, G, Com.npass, Com.nscale, nO) ;
            vs_DenseMult (z, b, V, n, nO) ;
            vs_free (V) ;
            /* (n, nO) = number of (rows, cols) in orthogonal matrix V
               z  = V*b = eigenvector estimate */
        }
 
        /* compute square of length of z (eigenvector) */
        s = vs_dot (z, z, n) ;
        /* r_e = sqrt(RE_SCALAR*n), s = r_e/||z||  */
        s = sqrt(RE_SCALAR*n/s) ;
        /* multiply z by s to obtain a vector whose length is r_e */
        vs_scale (z, s, n) ;
        /* translate z by the centroid c = 1/2 of the box constraint */
        vs_translate (z, VSHALF, n) ;
 
        vs_free (d) ;
        vs_free (u) ;
        vs_free (A) ;
        vs_free (Com.wi) ;
        vs_free (G->gs) ;
        vs_free (G->gi) ;
        vs_free (G->gj) ;
        vs_free (G->gk) ;
        vs_free (G->gf) ;
        vs_free (G->gx) ;
        vs_free (G->gy) ;
        vs_free (G->gz) ;
    }
    else /* VSP-based start guess */
    {
        for (i = 0; i < n; i++)
        {
            if ( x [i] == 2 ) z [i] = VSHALF ;
            else              z [i] = (VSFLOAT) x [i] ;
        }
    }
    /* Perform ngp gradient projection steps where the gradient is
       normalized to have length r_n = sqrt(n)/4. The projection is
       computed using napheap.  Here is the setup for napheap.
       Our version of napheap sets the lower and upper bounds to 0 and 1,
       so these parameters are constants, not variables.  */
    NAPdata *napdata ;
    VSFLOAT sqrtn, inv_n, *X ;
    sqrtn = 1.25*sqrt((VSFLOAT) n) ;
    inv_n = VSONE/n ;
    napdata = napheap_setup () ;
    napdata->c = z ; /* z is projected onto feasible set for (GB) */
    napdata->n = n ;
    if ( PrintLevel )
    {
        printf ("n = %i\n", n) ;
        printf ("m = %i\n", (int) m) ;
    }
    /* for (i = 0; i < n; i++)
    {
        printf ("a[%i] = %i ;\n", i, (int) a [i]) ;
    }*/
    napdata->a = a ;
    napdata->blo = m ;
    napdata->bhi = m ;
    napdata->x = X = (VSFLOAT *) vs_malloc (&status,  n, sizeof (VSFLOAT)) ;
    if ( status ) return (status) ;
    int const ngp = Parm->ngp ;

    /* evaluate the linear term = A*(1/2) in the objective of (GB). */
    A12 = (VSFLOAT *) vs_malloc (&status, n, sizeof (VSFLOAT)) ;
    if ( status ) return (status) ;
    vs_A1mult (A12, Ap, Ax, n, E_is_one) ;
    /* optionally project starting point z on feasible set for bisection */
    if ( Parm->initial_project )
    {
        status = napheap (napdata) ;
        /* move projection from X back to z */
        vs_copyx (z, X, n) ;
    }

    /* perform ngp gradient projection steps */
    for (i = 0; i < ngp; i++)
    {
        /* compute b = -A*z */
        vs_Amult (b, z, Ap, Ai, Ax, n) ;
        /* compute b += A12 ; */
        vs_vecadd (b, b, A12, n) ;

        /* determine scaling factor for the gradient */
        s = vs_dot (b, b, n) ;
        s = sqrt(0.125*m/s) ; /* sqrt((n/16)/(b'*b)) */
#if 0
        /* Decreasing the step size based on the distance to the surface of
           the sphere with center c and radius sqrt(n)/2 did not work as
           well. The larger fixed stepsize given above was better. */
        else /* i > 0 */
        {
             /* use 1.25*(sqrt(n)/2 - ||z-c||) since ||z-c|| approx
                         (sqrt(n)/2 at optimality */
             t = 0 ;
             for (i = 0; i < n; i++)
             {
                 s = z [i] - VSHALF ;
                 t += s*s ;
             }
             t = sqrtn*(0.5 - sqrt(inv_n*t)) ; /* desired length of grad vec */
             s = vs_dot (b, b, n) ;
             s *= t*sqrt(VSONE/s) ;
        }
#endif
        /* multiply b by -s to obtain a vector along negative gradient
           whose length is sqrt(n)/4 */
        vs_scale (b, -s, n) ;
        /* add the negative gradient vector to the current iterate */
        vs_vecadd (z, z, b, n) ;
        /*for (i = 0; i < n; i++)
        {
            printf ("y[%i] = %25.16e ;\n", i, z [i]) ;
        }
        fflush(stdout) ; */

        /* project z onto feasible set for (GB) */
        status = napheap (napdata) ;

        if ( PrintLevel )
        {
            VSFLOAT err ;
            err = napheap_check (&r, &s, &t, napdata) ;
            printf ("err: %e errb: %e errB: %e errc: %e\n", err, r, s, t) ;
            if ( err > 1.e-6 )
            {
                printf ("NAPHEAP error: %e\n", err) ;
                exit(1) ;
            }
        }

        /* move projection from X back to z */
        vs_copyx (z, X, n) ;
    }
    napheap_terminate (&napdata) ;
    
    /* Convert a feasible point for (GB) into feasible point on the boundary
       of the feasible set, where z [i] = 0 or 1 for n-1 choices of i in
      [0, n). The objective value for the boundary point <= objective value
      at the initial feasible point. */
    vs_boundary (z, Ai, Ap, Ax, a, n, A12, b, PrintLevel) ;

    if ( PrintLevel ) n0 = n1 = ns = cutedge = 0 ;
    p = 0 ;
    for (j = 0; j < n; j++)
    {
        if ( z [j] == 0 )
        {
            x [j] = 0 ;
            if ( PrintLevel )
            {
                if ( vw_is_one ) n0++ ;
                else             n0 += vw [j] ;
                VSLINT const l = Ap [j+1] ;
                for (; p < l; p++)
                {
                    if ( Ai [p] == j ) continue ;
                    if ( z [Ai [p]] == 1 ) cutedge++ ;
                }
            }
        }
        else if ( z [j] == 1 )
        {
            x [j] = 1 ;
            if ( PrintLevel )
            {
                if ( vw_is_one ) n1++ ;
                else             n1 += vw [j] ;
                /* look for z [j] = 0 */
                VSLINT const l = Ap [j+1] ;
                for (; p < l; p++)
                {
                    if ( Ai [p] == j ) continue ;
                    if ( z [Ai [p]] == 0 ) cutedge++ ;
                }
            }
        }
        else
        {
            if ( PrintLevel )
            {
                if ( vw_is_one ) ns++ ;
                else             ns += vw [j] ;
            }
            x [j] = 2 ;
        }
    }
    if ( PrintLevel )
    {
        cutedge /= 2 ;
        printf ("n0: %i n1: %i ns: %i cutedge: %i\n", n0, n1, ns, cutedge) ;
    }
    
    vs_separator
    (
        x, Ei, Ep, n, vw,
        Parm->try_both_covers, 0,
        NULL, NULL, NULL, NULL, NULL
    ) ;

    if ( PrintLevel )
    {
        n0 = n1 = ns = 0 ;
        for (j = 0; j < n; j++)
        {
            if      ( x [j] == 0 )
            {
                if ( vw_is_one ) n0++ ;
                else             n0 += vw [j] ;
            }
            else if ( x [j] == 1 )
            {
                if ( vw_is_one ) n1++ ;
                else             n1 += vw [j] ;
            }
            else
            {
                if ( vw_is_one ) ns++ ;
                else             ns += vw [j] ;
            }
        }
        printf ("n0: %i n1: %i ns: %i\n", n0, n1, ns) ;
    }

    if ( a != NULL ) vs_free (a) ;
    if ( Ax != NULL ) vs_free (Ax) ;
    vs_free (b) ;
    vs_free (z) ;
    vs_free (X) ;
    vs_free (Ai) ;
    vs_free (Ap) ;
    vs_free (A12) ;
    return (status) ;
}

/* =========================================================================
   ============================== vs_malloc ================================
   ========================================================================= */
void * vs_malloc
(
    int *status,
    VSINT     n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc ((int64_t) n * size) ;
#else
        p = malloc ((int64_t) n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                   *status = VS_OUT_OF_MEMORY ;
    else if ( *status != VS_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

/* =========================================================================
   ============================== vs_lmalloc ===============================
   ========================================================================= */
void * vs_lmalloc
(
    int *status,
    VSLINT     n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc (n * size) ;
#else
        p = malloc (n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                   *status = VS_OUT_OF_MEMORY ;
    else if ( *status != VS_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

/* =========================================================================
   ============================== vs_Lmalloc ===============================
   ========================================================================= */
void * vs_Lmalloc
(
    int *status,
    int64_t   n,
    int    size
)
{
    void *p ;
    if ( n > 0 )
    {
#ifdef MATLAB_MEX_FILE
        p = mxMalloc (n * size) ;
#else
        p = malloc (n * size) ;
#endif
    }
    else
    {
        return (NULL) ;
    }
    if      ( p == NULL )                   *status = VS_OUT_OF_MEMORY ;
    else if ( *status != VS_OUT_OF_MEMORY ) *status = 0 ;
    return (p) ;
}

/* =========================================================================
   ================================ vs_free ================================
   ========================================================================= */
void vs_free
(       
    void * p 
)       
{           
            
#ifdef MATLAB_MEX_FILE
    mxFree (p) ;
#else
    free (p) ;
#endif  
            
}

void vs_printAMATLAB
(
    VSINT   const ncol, /* number of cols in A */
    VSLINT  const  *Ap, /* size ncol+1, column pointers */
    VSINT   const  *Ai, /* size Ap [ncol], row indices for A */
    VSFLOAT const  *Ax, /* size Ap [ncol], numerical entries of A */
    char           *what  /* name of the matrix */
)
{
    VSINT j ;
    VSLINT p, q ;
    printf ("Temp = [\n") ;
    for (j = 0; j < ncol; j++)
    {
        q = Ap [j+1] ;
        {
            for (p = Ap [j]; p < q; p++)
            {
                if ( Ax == NULL )
                {
                    printf ("%i %i 1\n", Ai [p]+1, j+1) ;
                }
                else
                {
                    printf ("%i %i %25.15e\n",
                        Ai [p]+1, j+1, Ax [p]) ;
                }
            }
        }
    }
    printf ("] ;\n") ;
    printf ("%s = sparse (Temp(:, 1), Temp(:, 2), Temp(:, 3)) ;\n", what) ;
}

/* =========================================================================
   === vs_Amult ==========================================================
   =========================================================================
   Compute A * v and store it in Av
   ========================================================================= */
void vs_Amult
(
    VSFLOAT     *Av, // Stores the result, A * v
    VSFLOAT      *v, // Vector being multiplied by A
    VSLINT      *Ap, // Row pointers
    VSINT       *Ai, // Column indices
    VSFLOAT     *Ax, // Matrix entries, note that Ax < 0
    VSINT         n  // Number of columns
)
{
#ifdef USE_MKL
    VSFLOAT alpha, beta ;
    alpha = 1 ;
    beta  = 0 ;

    // Wrap the arrays in an mkl struct. Lightweight, no data copying
    sparse_matrix_t A ;
    struct matrix_descr descr = { SPARSE_MATRIX_TYPE_GENERAL, 0, 0 } ;

    mkl_sparse_d_create_csr
        (&A, SPARSE_INDEX_BASE_ZERO, n, n, Ap, Ap + 1, Ai, Ax) ;

    // Perform the update y = alpha * Av + beta * y
    mkl_sparse_d_mv
        (SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, v, beta, Av) ;
    
    // Cleanup the mkl struct
    // This does NOT destroy Ai, Ap, Ax arrays which were passed to it
    mkl_sparse_destroy(A) ;

#else

    VSINT j, j0 ;
    VSLINT p = 0 ;
    VSFLOAT t ;

    if ( Ax == NULL )
    {
        j0 = 0 ;
        for (j = 1 ; j <= n ;)
        {
            VSLINT const l = Ap [j] ;
            t = VSZERO ;
            for (; p < l; p++)
            {
                t -= v[Ai[p]] ;
            }
            Av[j0] = t ;
            j0 = j++ ;
        }
    }
    else
    {
        j0 = 0 ;
        for (j = 1 ; j <= n ;)
        {
            VSLINT const l = Ap[j] ;
            t = 0 ;
            for (; p < l; p++)
            {
                t  += Ax[p] * v[Ai[p]] ; /* Ax < 0 */
            }
            Av[j0] = t ;
            j0 = j++ ;
        }
    }
#endif
}

/* =========================================================================
   === vs_A1mult ===========================================================
   =========================================================================
   Compute A * 1/2 and store it in A12 (note: original A not A bar = -A)
   ========================================================================= */
void vs_A1mult
(
    VSFLOAT   *A12, // Stores the product
    VSLINT     *Ap, // Row pointers
    VSFLOAT    *Ax, // Matrix entries, Ax < 0, insert minus sign if E_is_one
    VSINT        n, // Number of columns
    int   E_is_one  // Flag (set to 1 for binary matrix, 0 otherwise) 
)
{
    VSINT j, j0 ;
    VSFLOAT t ;
    VSLINT p = 0 ;

    j0 = 0 ;
    if ( E_is_one )
    {
        for (j = 1 ; j <= n ;)
        {
            A12 [j0] = VSHALF*(Ap [j] - p) ;
            p = Ap [j] ;
            j0 = j++ ;
        }
    }
    else
    {
        for (j = 1 ; j <= n ;)
        {
            VSLINT const l = Ap[j] ;
            t = 0 ;
            for (; p < l; p++)
            {
                t -= Ax[p] ; /* minus sign needed to recover original A */
            }
            A12 [j0] = VSHALF*t ;
            j0 = j++ ;
        }
    }
}

/* =========================================================================
   === vs_sort_cols ======================================================
   =========================================================================
   Given a symmetric matrix, output the same matrix but with sorted columns
   ========================================================================= */
int vs_sort_cols
(
    VSLINT      *Ap, /* given column pointers for symmetric matrix */
    VSINT       *Ai, /* sorted row indices (spaced malloc'd by user) */
    VSINT       *Ax, /* sorted numerical entries (spaced malloc'd by user) */
    VSINT const *Bi, /* given row indices (unsorted) */
    VSINT const *Bx, /* given entries (unsorted), Bx = NULL => Ax = NULL */
    VSINT         n  /* number of rows in matrix */
)
{
    VSINT i, j ;
    VSLINT p, *Bp ;
    Bp = malloc ((n+1)*sizeof (VSLINT)) ;
    if ( Bp == NULL ) return (1) ;
    int const AxExists = (Bx != NULL) ? 1 : 0 ;
    for (j = 0; j < n; j++) Bp [j] = Ap [j] ;
    p = 0 ;
    for (j = 0; j < n; j++)
    {
        VSLINT const q = Ap [j+1] ;
        for (; p < q; p++)
        {
            VSLINT const pp = Bp [Bi [p]]++ ;
            Ai [pp] = j ;
            if ( AxExists ) Ax [pp] = Bx [p] ;
        }
    }
    free (Bp) ;
    return (0) ;
}

/* =========================================================================
   === vs_daxpy ==========================================================
   =========================================================================
   Evaluate y = a*x + y
   ========================================================================= */
void vs_daxpy
(
    VSFLOAT *y,
    VSFLOAT const *x,
    VSFLOAT const  a,
    VSINT   const  n
)
{
#ifdef USE_MKL
    // The 1's are the x and y increments
    cblas_daxpy(n, a, x, 1, y, 1) ;
#else
    VSINT i ;
    for(i = 0 ; i < n ; i++) y[i] += a*x[i] ;
#endif
}

VSFLOAT vs_max
(
    VSFLOAT const *x, /* VSFLOAT vector */
    VSINT   const  n  /* VSINT length of vector */
)
{
    VSINT i, n5, t ;
    t = 0 ;
    n5 = n % 5 ;

    for (i = 0; i < n5; i++) if ( t < fabs(x [i]) ) t = fabs(x [i]) ;
    for (; i < n; i += 5)
    {
        if ( t < x [i]   ) t = fabs(x [i])   ;
        if ( t < x [i+1] ) t = fabs(x [i+1]) ;
        if ( t < x [i+2] ) t = fabs(x [i+2]) ;
        if ( t < x [i+3] ) t = fabs(x [i+3]) ;
        if ( t < x [i+4] ) t = fabs(x [i+4]) ;
    }
    return (t) ;
}

/* =========================================================================
   === vs_norm ===========================================================
   =========================================================================
   Euclidean norm of a vector
   ========================================================================= */
VSFLOAT vs_norm
(
    VSFLOAT const *x, /* return Euclidean norm of x */
    VSINT   const  n  /* dimension */
)
{
    return (sqrt (vs_dot (x, x, n))) ;
}

void vs_initx
(
    VSFLOAT      *x,  /* array to be initialized */
    VSFLOAT const s,  /* scalar */
    VSINT  const n   /* length of x */
)
{
    VSINT j, n5 ;
    VSFLOAT *xj ;
    if ( n == 0 ) return ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

void vs_initi
(
    VSINT      *x,  /* array to be initialized */
    VSINT const s,  /* scalar */
    VSINT const n   /* length of x */
)
{
    VSINT j, n5 ;
    VSINT *xj ;
    if ( n == 0 ) return ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

void vs_initLx
(
    VSFLOAT      *x,  /* array to be initialized */
    VSFLOAT const s,  /* scalar */
    VSLINT  const n   /* length of x */
)
{
    VSINT j, n5 ;
    VSFLOAT *xj ;
    if ( n == 0 ) return ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] = s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
        *(xj++) = s ;
    }
}

/* =========================================================================
   === vs_copyx ============================================================
   =========================================================================
   Copy VSFLOAT vector y into vector x
   ========================================================================= */
void vs_copyx
(
    VSFLOAT       *x, /* output of copy */
    VSFLOAT const *y, /* input of copy */
    VSINT   const  n  /* length of vectors */
)
{
    if ( (y == x) || (y == NULL) ) return ;
    {
        VSINT i, n5 ;
        n5 = n % 5 ;
        for (i = 0; i < n5; i++) x [i] = y [i] ;
        for (; i < n; i += 5)
        {
            x [i]   = y [i] ;
            x [i+1] = y [i+1] ;
            x [i+2] = y [i+2] ;
            x [i+3] = y [i+3] ;
            x [i+4] = y [i+4] ;
        }
    }
}

void vs_scale
(
    VSFLOAT       *x,  /* scaled array */
    VSFLOAT const  s,  /* scale */
    VSINT   const  n   /* length of x */
)
{
    VSFLOAT *xj ;
    if ( s == 1 ) return ; /* nothing to do */
    VSINT j, n5 ;
    n5 = n % 5 ;
    for (j = 0; j < n5; j++) x [j] *= s ;
    xj = x+j ;
    for (; j < n; j += 5)
    {
        *(xj++) *= s ;
        *(xj++) *= s ;
        *(xj++) *= s ;
        *(xj++) *= s ;
        *(xj++) *= s ;
    }
}

/* =========================================================================
   === vs_project ==========================================================
   =========================================================================
   Replace x by its projection onto the space orthogonal to a
   ========================================================================= */
void vs_project
(
    VSFLOAT              *x, /* projection */
    VSFLOAT const        *a, /* project onto orthogonal complement of a */
    VSFLOAT const        aa, /* 1/a'a */
    int     const vw_is_one,
    VSINT   const         n  /* dimension */
)
{
    VSINT j, n5 ;
    VSFLOAT t, ax ;
    n5 = n%5 ;
    t = VSZERO ;
    if ( vw_is_one )
    {
        for (j = 0; j < n5; j++)
        {
            t += x [j] ;
        }
        for (; j < n; j += 5)
        {
            t += x [j] + x[j+1] + x[j+2] + x[j+3] + x[j+4] ;
        }
        t *= aa ;
        for (j = 0; j < n5; j++)
        {
            x [j] -= t ;
        }
        for (; j < n;)
        {
            x [j++] -= t ;
            x [j++] -= t ;
            x [j++] -= t ;
            x [j++] -= t ;
            x [j++] -= t ;
        }
    }
    else /* a != 1 */
    {
         ax = vs_dot (a, x, n)*aa ;
         vs_daxpy (x, a, -ax, n) ;
    }
}

/* =========================================================================
   === vs_dot ==============================================================
   =========================================================================
   Return the dot product a'*x
   ========================================================================= */
VSFLOAT vs_dot
(
    VSFLOAT const *a,
    VSFLOAT const *x,
    VSINT   const  n
)
{
    VSINT n5, j ;
    VSFLOAT t ;
    t = VSZERO ;
    n5 = n%5 ;
    if ( a == NULL )
    {
        for (j = 0; j < n5; j++) t += x [j] ;
        for (; j < n; j += 5)
        {
            t += x [j] + x [j+1] + x [j+2] + x [j+3] + x [j+4] ;
        }
    }
    else
    {
        for (j = 0; j < n5; j++) t += a [j]*x [j] ;
        for (; j < n; j += 5)
        {
            t += a[j]*x[j] + a[j+1]*x[j+1] + a[j+2]*x[j+2] + a[j+3]*x[j+3]
                           + a[j+4]*x[j+4] ;
        }
    }
    return (t) ;
}

/* =========================================================================
   === vs_translate ========================================================
   =========================================================================
   Set new += t
   ========================================================================= */
void vs_translate
(
    VSFLOAT      *x,
    VSFLOAT const t,
    VSINT   const n
)
{
    VSINT j ;
    VSINT const n5 = n%5 ;
    for (j = 0; j < n5; j++)
    {
        x [j] += t ;
    }

    VSFLOAT *xj ;
    xj = x+j ;

    for (; j < n; j += 5)
    {
        *(xj++) += t ;
        *(xj++) += t ;
        *(xj++) += t ;
        *(xj++) += t ;
        *(xj++) += t ;
    }
}

/* ==========================================================================
    === vs_vecadd ===========================================================
    Set x = y + b 
============================================================================= */
void vs_vecadd
(
     VSFLOAT       *x,
     VSFLOAT       *y,
     VSFLOAT const *b,
     VSINT   const  n
)
{
     VSINT j ;
     VSINT const n5 = n%5 ;
     if ( x == y )
     {
         for (j = 0; j < n5; j++)
         {
             x [j] += b [j] ;
         }

         for (; j < n; j += 5)
         {
             x [j]   += b [j] ;
             x [j+1] += b [j+1] ;
             x [j+2] += b [j+2] ;
             x [j+3] += b [j+3] ;
             x [j+4] += b [j+4] ;
         }
     }
     else
     {
         for (j = 0; j < n5; j++)
         {
             x [j] = y [j] + b [j] ;
         }

         for (; j < n; j += 5)
         {
             x [j]   = y [j]   + b [j] ;
             x [j+1] = y [j+1] + b [j+1] ;
             x [j+2] = y [j+2] + b [j+2] ;
             x [j+3] = y [j+3] + b [j+3] ;
             x [j+4] = y [j+4] + b [j+4] ;
         }
     }
}

/* =========================================================================
   ============================= vs_boundary ===============================
   =========================================================================
       ________________________________________________________________
      |    Consider the following version of the graph bisection       |
      |    problem:                                                    |
      |                                                                |
      |                           1                                    |
      |    (GB)       min f(x) := - (1-x)'Ax                           |
      |                           2                                    |
      |               subject to a'x = m, 0 <= x <= 1                  |
      |                                                                |
      |    n = number of vertices in an undirected graph, a > 0,       |
      |    m = (sum a_i)/2, A = D + E where where E is the edge        |
      |    weighted adjacency matrix for graph (diagonal of E is zero),|
      |    and D is a diagonal matrix with d_{jj} the maximum element  |
      |    of E in column j. It is assumed that A has no zero columns. |
      |    Given a feasible point x, the code outputs a new feasible   |
      |    point x* with f(x*) <= f(x), and with x [i] = 0 or 1        |
      |    for n-1 distinct i in [0, n). The final value of i should   |
      |    satisfy 0 <= x [i] <= 1, however, both inequalities are     | 
      |    potentially strict.                                         |              |                                                                |
      |    Algorithm:                                                  |
      |                                                                |
      |    Search for two components of x with                         |
      |                                                                |
      |        0 < x [i] < 1, 0 < x [j] < 1, and i < j                 |
      |                                                                |
      |    Let p denote a vector with all components zero except for   |
      |    p [i] = a [j] and p [j] = -a [i]. For this choice of p,     |
      |    we have a'p = 0. Hence, x + s * p satisfies the linear      |
      |    constraint of (GB) for any choice of s. By a Taylor         |
      |    expansion of f around x, we have                            |
      |                                                                |
      |        f(x + s * p) = f(x) + s * grad f(x) p - 0.5 s*s p'Ap    |
      |                                                                |
      |    By the choice for D, it follows that p'Ap >= 0 where the    |
      |    inequality is strict if a [i] != a [j] or d [i] > a_{ij} or |
      |    d [j] > a_{ij}. We choose s as follows:                     |
      |                                                                |
      |        Case 1. If grad f (x) p > 0, then take s < 0, and if    |
      |                   grad f (x) p < 0, then take s > 0.           |
      |                Choosing s in this way, the linear term is      |
      |                strictly negative and the second derivative     |
      |                term is nonpositive. Take the largest |s|       |
      |                such that x + s * p is feasible for (GB).       |
      |                                                                |
      |        Case 2. If grad f (x) p = 0, then choose |s| as large as|
      |                possible with x + s * p feasible for (GB)       |
      |                                                                |
      |        In case 1, the linear term is strictly negative and the |
      |        quadratic term is nonnegative. We take the largest |s|  |
      |        for which x + s * p is feasible. In case 2, the linear  |
      |        term is zero and the quadratic term remains er indices  |
      |        greater than j for a new i and j. If the former holds,  |
      |        then either the old i or old j becomes the new i and    |
      |        we search over indices greater than j for a new j.      |
      |        Besides the update to i and j, the gradient grad f(x)   |
      |        should be updated after the move in x. The process ends |
      |        when a new j does not exist.                            |
      |                                                                |
      |                   Copyright by William W. Hager                |
      |________________________________________________________________|

       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|*/

void vs_boundary /* return original cost and its improvement */
(
    VSFLOAT           *x, /* feasible point for (GB), output described above */
    VSINT            *Ai, /* row numbers */
    VSLINT           *Ap, /* column pointers */
    VSFLOAT          *Ax, /* edge weights,   NULL => all entries in Ax are 1 */
    VSFLOAT           *a, /* vertex weights, NULL => all entries in a are 1 */
    VSINT              n, /* problem dimension */
    VSFLOAT         *A12, /* linear term in objective of (GB) */
    VSFLOAT           *g, /* work array for storing gradient vector */
    int const PrintLevel
)
{
    int dumpi, dumpj, two_hit ;
    VSFLOAT pi, pj, s1, s2, s3, s4, sposmax, snegmax,
            xi, xj, tlinear ;
    VSINT i, j ;
    VSLINT p ;

    int const E_is_one = (Ax == NULL) ? VSTRUE : VSFALSE ;
    /* evaluate gradient at current x */
    /* compute g = -A*x */
    vs_Amult (g, x, Ap, Ai, Ax, n) ;

    /* int status = 0 ;
    VSFLOAT *cost ;
    cost = vs_malloc (&status, 2, sizeof (VSFLOAT)) ;
    cost [0] = vs_dot (A12, x, n) + 0.5*vs_dot (x, g, n) ;
    VSFLOAT t, *gnew ;
    VSINT k ;
    gnew = malloc (n*sizeof(VSFLOAT)) ; */

    /* find index > -1 */
    i = vs_interior (-1, x, n) ;
    if ( i == VSEMPTY )
    {
        /* x has no interior components */
        /* cost [1] = VSZERO ;
        return (cost) ; */

        return ;
    }
    j = vs_interior (i, x, n) ;
    if ( j == VSEMPTY )
    {
        /* done since x has at most one fractional component */
        /* cost [1] = VSZERO ;
        return (cost) ;*/

        return ;
    }
    /* compute g += A12 ; */
    vs_vecadd (g, g, A12, n) ;
    int const a_is_one = (a == NULL) ? VSTRUE : VSFALSE ;
    two_hit = VSFALSE ;
    dumpi = VSFALSE ;
    dumpj = VSFALSE ;
    while ( 1 )
    {
        /* components i and j of x are fractional */
        VSFLOAT const ai = (a_is_one) ? VSONE : a [i] ;
        VSFLOAT const aj = (a_is_one) ? VSONE : a [j] ;
        VSFLOAT const xi = x [i] ;
        VSFLOAT const xj = x [j] ;
        /* pi = aj, pj = -ai */
        tlinear = g [i]*aj - g [j]*ai ;

        /* if tlinear = 0, then determine whether a positive or a negative s
           step is larger in absolute value, and take the larger step */
        if ( tlinear == VSZERO )
        {
            /* pos s step */
            s1 = (1-xi)/aj ; /* xi hits upper bound */
            s2 = xj/ai ;     /* xj hits lower bound */
            if ( s1 < s2 )
            {
                sposmax = s1 ;
            }
            else sposmax = s2 ;

            /* neg s step */
            s3 = xi/aj ;       /* xi hits lower bound */
            s4 = (1-xj)/ai ;   /* xj hits upper bound */
            if ( s3 < s4 ) snegmax = s3 ;
            else           snegmax = s4 ;

            if ( sposmax > snegmax )
            {
                if ( sposmax == s1 )
                {
                    /* s is positive and xi hits upper bound */
                    pi = VSONE - xi ;
                    x [i] = VSONE ;
                    pj = -s1*ai ;
                    x [j] += pj ;
                    if ( x [j] <= VSZERO )
                    {
                        if ( PrintLevel )
                        {
                            printf ("x [%i]: %e violates lower bound\n",
                                     j, x [j]) ;
                        }
                        two_hit = VSTRUE ;
                        x [j] = VSZERO ;
                    }
                    else dumpi = VSTRUE ;
                }
                else /* sposmax = s2 */
                {
                    /* s is positive and xj hits lower bound */
                    pj = -xj ;
                    x [j] = VSZERO ;
                    pi = s2*aj ;
                    x [i] += pi ;
                    if ( x [i] >= VSONE )
                    {
                        if ( PrintLevel )
                        {
                            printf ("x [%i]: %e violates upper bound\n",
                                     i, x [i]) ;
                        }
                        two_hit = VSTRUE ;
                        x [i] = VSONE ;
                    }
                    else dumpj = VSTRUE ;
                }
            }
            else
            {
                if ( snegmax == s3 )
                {
                    /* s is negative and xi hits lower bound */
                    pi = -xi ;
                    x [i] = VSZERO ;
                    pj = s3*ai ;
                    x [j] += pj ;
                    if ( x [j] >= VSONE )
                    {
                        if ( PrintLevel )
                        {
                            printf ("x [%i]: %e violates upper bound\n",
                                     j, x [j]) ;
                        }
                        two_hit = VSTRUE ;
                        x [j] = VSONE ;
                    }
                    else dumpi = VSTRUE ;
                }
                else /* s is negative and xj hits upper bound */
                {
                    pj = VSONE - xj ;
                    x [j] = VSONE ;
                    pi = -s4*aj ;
                    x [i] += pi ;
                    if ( x [i] <= VSZERO )
                    {
                        if ( PrintLevel )
                        {
                            printf ("x [%i]: %e violates lower bound\n",
                                     i, x [i]) ;
                        }
                        two_hit = VSTRUE ;
                        x [i] = VSZERO ;
                    }
                    else dumpj = VSTRUE ;
                }
            }
        }
        else if ( tlinear > VSZERO )
        {
            /* neg s step */
            s3 = xi/aj ;       /* xi hits lower bound */
            s4 = (1-xj)/ai ;   /* xj hits upper bound */
            if ( s3 < s4 )     /* xi hits lower bound s3 smaller step */
            {
                pi = -xi ;
                x [i] = VSZERO ;
                pj = s3*ai ;
                x [j] += pj ;
                if ( x [j] >= VSONE )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %e violates upper bound\n", j, x [j]) ;
                    }
                    two_hit = VSTRUE ;
                    x [j] = VSONE ;
                }
                else dumpi = VSTRUE ;
            }
            else               /* xj hits upper bound, s4 smaller step */
            {
                pj = VSONE - xj ;
                x [j] = VSONE ;
                pi = -s4*aj ;
                x [i] += pi ;
                if ( x [i] <= VSZERO )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %e violates lower bound\n", i, x [i]) ;
                    }
                    two_hit = VSTRUE ;
                    x [i] = VSZERO ;
                }
                else dumpj = VSTRUE ;
            }
        }
        else     /* tlinear < 0 */
        {
            /* pos s step */
            s1 = (1-xi)/aj ; /* xi hits upper bound */
            s2 = xj/ai ;     /* xj hits lower bound */
            if ( s1 < s2 )   /* xi hits upper bound, s1 is smaller step */
            {
                pi = VSONE - xi ;
                x [i] = VSONE ;
                pj = -s1*ai ;
                x [j] += pj ;
                if ( x [j] <= VSZERO )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %e violates lower bound\n", j, x [j]) ;
                    }
                    two_hit = VSTRUE ;
                    x [j] = VSZERO ;
                }
                else dumpi = VSTRUE ;
            }
            else            /* xj hits lower bound, s2 is smaller step */
            {
                pj = -xj ;
                x [j] = VSZERO ;
                pi = s2*aj ;
                x [i] += pi ;
                if ( x [i] >= VSONE )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %e violates upper bound\n", i, x [i]) ;
                    }
                    two_hit = VSTRUE ;
                    x [i] = VSONE ;
                }
                else dumpj = VSTRUE ;
            }
        }

        /* update gradient to account for change in xi */
        VSLINT const l1 = Ap [i+1] ;
        if ( E_is_one )
        {
            for (p = Ap [i]; p < l1; p++)
            {
                g [Ai [p]] -= pi ;
            }
        }
        else
        {

            for (p = Ap [i]; p < l1; p++)
            {
                g [Ai [p]] += pi*Ax [p] ;
            }
        }
        /* update gradient to account for change in xj */
        VSLINT const l2 = Ap [j+1] ;
        if ( E_is_one )
        {
            for (p = Ap [j]; p < l2; p++)
            {
                g [Ai [p]] -= pj ;
            }
        }
        else
        {
            for (p = Ap [j]; p < l2; p++)
            {
                g [Ai [p]] += pj*Ax [p] ;
            }
        }

        /* vs_Amult (gnew, x, Ap, Ai, Ax, n) ;
        vs_vecadd (gnew, gnew, A12, n) ;
        t = 0 ;
        for (k = 0; k < n; k++)
        {
            t += fabs(gnew [k] - g [k]) ;
        }   
        printf ("error in gradient: %e\n", t) ; */

        /* update i and j if two_hit is TRUE */
        if ( two_hit )
        {
            /* find first index > j */
            i = vs_interior (j, x, n) ;
            if ( i == VSEMPTY )
            {
                /* x has no interior components */
                break ;
            }
            // Buggy
            // j = vs_interior (j, x, n) ;
            j = vs_interior (i, x, n) ;
            if ( j == VSEMPTY )
            {
                /* done since x has at most one fractional component */
                break ;
            }
            two_hit = VSFALSE ;
        }
        else if ( dumpi )
        {
            i = j ;
            j = vs_interior (j, x, n) ;
            if ( j == VSEMPTY ) break ;
            dumpi = VSFALSE ;
        }
        else /* dumpj */
        {
            j = vs_interior (j, x, n) ;
            if ( j == VSEMPTY ) break ;
            dumpj = VSFALSE ;
        }
    }
    /* cost [1] = cost [0] - 0.5*(vs_dot (A12, x, n) + vs_dot (x, g, n)) ;
    if ( (cost [1]/cost [0]) < -1.e-12 )
    {
        printf ("boundary code hurts objective value\n") ;
        fflush(stdout) ;
        exit(1) ;
    }
    free (cost) ; */

    return ;
}

VSINT vs_interior
(
   VSINT    ilower, /* search x for indices > ilower and < n */
   VSFLOAT      *x, /* search for 0 < x [i] < 1 */
   VSINT         n  /* dimension of x */
)
{
    VSINT i ;
    for (i = ilower+1; i < n; i++)
    {
        if ( (x [i] > VSZERO) && (x [i] < VSONE) ) break ;
    }
    if ( i < n ) return (i) ;
    else         return (VSEMPTY) ;
}

/* ==========================================================================
   === vs_house =============================================================
   ==========================================================================
  
    Reduce A to tridiagonal form by Householder orthogonal similarity
    transformations:
                                 _
                                | d_1  u_1   0    0    .  .  .
                                |
             V'AV = T =         | u_1  d_2   u_2   0   .  .  .
                                |
                                |  0   u_2   d_3   0   .  .  .
                                |
                                |  .    .     .    .   .  .  .

   where V is computed as a product of Householder matrices. A is
   stored as a dense symmetric matrix, and V is stored as a dense matrix.
   The output matrix V overwrites the input matrix A.

n = size (A, 1) ;
d = zeros(n,1) ;
u = zeros(n,1) ;

if ( n == 1 ) % {

    d(1) = A(1,1) ;
    V(1) = 1 ;

elseif ( n == 2 ) % } {

    d(1) = A(1,1) ;
    d(2) = A(2,2) ;
    u(1) = A(2,1) ;
    V = zeros(2) ;
    V(1,1) = 1 ;
    V(2,2) = 2 ;

else % } {

    W = zeros(n,n) ; Householder vectors (stored in A)
    for j = 1:n-2 % {

        %-----------------------------------------------------------------------
        % compute v
        %-----------------------------------------------------------------------

        % accesses only the strictly lower triangular part of A:
        v = house (A(:,j),j+1) ;
        % v (1:j) is zero

        %-----------------------------------------------------------------------
        % compute column j+1 of W
        %-----------------------------------------------------------------------

        % column 1 of W is zero.  W is strictly lower triangular
        W(:,j+1) = v ;

        %-----------------------------------------------------------------------
        % compute x
        %-----------------------------------------------------------------------

        % since v (j) is zero, we can skip column j of A
        % also note that A (1:(j-1),(j+1):n) is zero
        % accesses just lower part of A to compute this
        x = zeros (n,1) ;
        x (j:n) = A (j:n,(j+1):n) * v ((j+1):n) ;

        %-----------------------------------------------------------------------
        % compute a
        %-----------------------------------------------------------------------

        % x (1:(j-1)) is zero and v (1:j) is zero
        a = .5 * ( (v ((j+1):n))' * x ((j+1):n)) ;

        %-----------------------------------------------------------------------
        % compute b
        %-----------------------------------------------------------------------

        b = zeros (n,1) ;
        b (j:n) = a * v (j:n) - x (j:n) ;

        %-----------------------------------------------------------------------
        % update A
        %-----------------------------------------------------------------------

        % just update the lower triangular part of A
        A(j:n,j:n) = A(j:n,j:n) + v(j:n)*b(j:n)' + b(j:n)*v(j:n)' ;

        %-----------------------------------------------------------------------
        % save diagonal and off-diagonal
        %-----------------------------------------------------------------------

        d(j) = A(j,j) ;
        u(j) = A(j+1,j) ;

    end % }
    d(n-1) = A(n-1,n-1) ;
    u(n-1) = A(n,n-1) ;
    d(n) = A(n,n) ;

    %---------------------------------------------------------------------------
    % project W to V
    %---------------------------------------------------------------------------

    V = zeros(n,n) ;
    for j = 1:n % {
        V(:,j) = qj (W, j) ;
    end % }

end % }

   ========================================================================== */

void vs_house
(
    VSFLOAT *A, /* input:  n-by-n dense matrix
                    output: tridiagonalization plus Household product in A */
    VSFLOAT *d, /* n-by-1 vector (output) */
    VSFLOAT *u, /* n-by-1 vector (output) */
    VSINT    n, /* matrix dimension */
    VScom  *com  /* common data */
)
{
    VSFLOAT a, hj, y, s, t, *Ak, *Aj, *x ;
    VSINT i, j, jp1, k ;

    if (n == 1)
    {

        d [0] = A [0] ;
        u [0] = 0 ;
        A [0] = 1 ;
        return ;

    }
    else if (n == 2)
    {

        d [0] = A [0] ;
        d [1] = A [3] ;
        u [0] = A [1] ;
        A [0] = 1 ;
        A [1] = 0 ;
        A [2] = 0 ;
        A [3] = 1 ;
        return ;

    }

    x = com->wx ;

     /* start Householder reduction for matrix of dimension >= 3 */
     /* initialize current column of A */
    Aj = A ;
    for (j = 0 ; j < n-2 ; j++)
    {
        d [j] = Aj [j] ;
        s = VSZERO ;
        jp1 = j+1 ;
        for (i = jp1; i < n; i++) s += Aj [i]*Aj [i] ;
        if ( s == VSZERO )
        {
            u [j] = VSZERO ;
            Aj += n ;
            continue ;   /* the Householder vector is 0 as is subcolumn of A */
        }

        /* scaling factor for Household matrix */
        s = sqrt (s) ;
        hj = Aj [jp1] ;
        t = 1/sqrt (s*(s + fabs (hj))) ;
        if ( hj >= 0 )
        {
            hj += s ;    /* hj is (j+1)st element of Household matrix */
            y = -s ;     /* y = u_j = (j+1)st element of (Householder * Aj) */
        }
        else
        {
            hj -= s ;
            y = s ;
        }
        u [j] = y ;

        /* ----------------------------------------------------------- */
        /* overwrite column j of A with Householder vector */
        /* ----------------------------------------------------------- */

        Aj [jp1] = hj ;
        for (i = jp1; i < n; i++) Aj [i] *= t ;
    
        /* x ((j+1):n) = A ((j+1):n,(j+1):n) * house ((j+1):n) */
        for (i = jp1; i < n; i++) x [i] = VSZERO ;
        Ak = Aj ;
        for (k = jp1 ; k < n ; k++)
        {
    	    /* Ak is column k of A */
            Ak += n ;
            /* diagonal */
            {
                VSFLOAT vk = Aj [k] ;
                VSFLOAT xk = Ak [k] * vk ;
    
                /* dot product with row k of A and saxpy with column k of A */
                for (i = k+1 ; i < n ; i++)
                {
                    VSFLOAT aki = Ak [i] ;
                    xk += aki * Aj [i] ;
                    x [i] += aki * vk ;
                }
                x [k] += xk ;
            }
        }

        /* -------------------------------------------------------------- */
        /* a = .5 * ( (v ((j+1):n))' * x ((j+1):n)) ; */
        /* -------------------------------------------------------------- */

        a = 0 ;
        for (i = jp1; i < n; i++) a += Aj [i] * x [i] ;
        a /= 2 ;

        /* ----------------------------------------------------- */
        /* b (j:n) = a * v (j:n) - x (j:n) ; write b back into x */
        /* ----------------------------------------------------- */

        for (i = jp1; i < n; i++)
        {
            x [i] = a*Aj [i] - x [i] ;
        }

        /* ----------------------------------------------------------------- */
        /* A(j+1:n,j+1:n) = A(j+1:n,j+1:n)
                            + v(j+1:n)*b(j+1:n)' + b(j+1:n)*v(j+1:n)' */
        /* ----------------------------------------------------------------- */

        /* note that x holds b */
        /* update only the lower triangular part of A */
        Ak = Aj ;
        for (k = jp1 ; k < n ; k++)
        {
            Ak += n ;
            {
                VSFLOAT vk = Aj [k] ;
                VSFLOAT xk = x [k] ;
                Ak [k] += 2.*(vk * xk) ;
                for (i = k+1; i < n; i++)
                {
                    Ak [i] += (Aj [i] * xk) + (x [i] * vk) ;
                }
            }
        }
        Aj += n ;
    }

    /* ------------------------------------------------------------------ */
    /* d(n-1) = A(n-1,n-1) ; u(n-1) = A(n,n-1) ; u(n) = 0 ; d(n) = A(n,n) */
    /* ------------------------------------------------------------------ */

    d [n-2] = Aj [n-2] ;
    u [n-2] = Aj [n-1] ;
    Aj += n ;           /* Aj = last column of matrix */
    d [n-1] = Aj [n-1] ;
    u [n-1] = 0 ;

    /* ------------------------------------------------------------------ */
    /* compute columns of V and store them in A */
    /* start with last column and work to first column */
    /* ------------------------------------------------------------------ */

    /* last column is special */
    Ak = Aj-(n+n) ;    /* first column with a Householder vector */
    t = Ak [n-1] ;
    Aj [n-2] =   - t*Ak [n-2] ;
    Aj [n-1] = 1 - t*Ak [n-1] ;
    for (j = n-4; j >= 0; j--)
    {
        /* compute the new element in Aj */
        Ak -= n ;
        t = VSZERO ;
        for (i = j+2; i < n; i++)
        {
            t += Aj [i]*Ak [i] ;
        }
        Aj [j+1] = -Ak [j+1]*t ;

        /* update the remaining elements in Aj */
        for (i = j+2; i < n; i++)
        {
            Aj [i] -= Ak [i]*t ;
        }
    }
    Aj [0] = 0 ;

    /* compute next to last through 2nd column of V */
    for (j = n-2; j > 0; j--)
    {
        Aj -= n ;
        /* startup, column of identity times preceding Householder */
        Ak = Aj - n ;
        t = Ak [j] ;
        for (i = j; i < n; i++) Aj [i] = -Ak [i]*t ;
        Aj [j]++ ;
        for (k = j-2; k >= 0; k--)
        {
            /* compute the new element in Aj */
            Ak -= n ;
            t = VSZERO ;
            for (i = k+2; i < n; i++)
            {
                t += Aj [i]*Ak [i] ;
            }
            Aj [k+1] = -Ak [k+1]*t ;

            /* update the remaining elements in Aj */
            for (i = k+2; i < n; i++)
            {
                Aj [i] -= Ak [i]*t ;
            }
        }
        Aj [0] = VSZERO ;
    }
    /* store the first columns of V = first column of identity */
    Aj -= n ;
    vs_initx (Aj, VSZERO, n) ;
    Aj [0] = 1 ;
}

/* ========================================================================== */
/* ==== vs_diag ============================================================= */
/* ========================================================================== */
/*
        Input: An n-dimensional tridiagonal matrix stored in arrays d and u

       Output: The eigenvalues along with the Givens rotations needed to
               diagonalize the matrix

    Algorithm: Inverse shifted QR algorithm with fast Givens rotations.
               The number of Givens rotations needed for the diagonalization
               is related to the number of QR iterations. This algorithm
               is an explicit implicit version of the QR method. It is
               explicit since the rotations are computed as in a
               QR factorization. It is implicit since each rotation
               is applied to both the right and left side of the matrix
               simultaneously, as is done with an implicit implementation
               of the QR method. An advantage of this explicit implicit
               algorithm is that the loss of shift phenomenon associated
               with the implicit QR algorithm does not occur since the
               rotations are computed from the explicit QR factorization.
               The speed associated with the implicit QR method is
               retained since each iteration entails a single pass over
               the matrix from upper left to lower right corner. In
               the fast Givens implementation, the matrix is stored as
               a product DTD where D is diagonal and T is tridiagonal.
               In the QR algorithm, the Givens rotations are applied to
               DTD - sigma I = D(T - sigma/D^2)D. In the explicit QR
               algorithm, the partially factored matrix has the form:

                         x  x  f  0  0  0  0
                         0  x  x  f  0  0  0
                         0  0  p  x  0  0  0
                         0  0  q  d  u  0  0
                         0  0  0  u  d  u  0
                         0  0  0  0  u  d  u

               Here f denotes fill, d and u are the original diagonal,
               subdiagonal, and superdiagonal elements, and p and q
               are the elements on the diagonal and subdiagonal which
               are used to generate the next Given rotation G. The
               rotation is designed to annihilate q while replacing
               p by sqrt(p*p + q*q). In performing the left and right
               multiplications by G, one needs to recall the structure
               of the matrix that is operated on:

                         d  u  0  0  0
                         u  d  u  a  0
                         0  u  d  u  0
                         0  a  u  d  u
                         0  0  0  u  d

               The matrix is tridiagonal with a bulge corresponding to
               the element a. In the code below, the iteration starts
               below the comment "perform one pass". The variables p, q,
               and a inside the subsequent code are identical to the
               p, q, and a above.

               An advantage of the fast Givens rotation is that several
               multiplications associated with the usual Givens rotation
               are eliminated as well as the square root. For details
               concerning fast Givens, see W. W. Hager, Applied Numerical
               Linear Algebra, pages 221-223. In the code that follows, the
               squares of the diagonal of the diagonal scaling matrix D are
               stored in the array sc (scaling array). As the algorithm
               progresses, the elements of sc tend to 0 since they are being
               multiplied by factors between 1/2 and 1.  When sc[i] drops
               below qr_lower, we multiply row i of T by sc[i] and reset
               sc[i] = 1. The code below rectifies an unflow problem
               which could occur in the analogous NAPACK code tdg.f.
               In tdg.f, we check sc at the start of each qr iteration.
               However, for a matrix with 5000 or more rows, an underflow
               could potentially occur within a QR iteration. Hence, in
               the new code, we check for small scales within the iteration.
               The code also checks to see if the matrix can be split due
               to a small off-diagonal element.
   ========================================================================== */
int vs_diag /* return: 0 (convergence tolerance satisfied)
                              -5 (algorithm did not converge)
                              -6 (insufficient space allocated for Givens) */
(
    VSFLOAT   *d, /* diagonal of tridiagonal matrix, output: eigenvalues */
    VSFLOAT   *u, /* superdiagonal of tridiagonal matrix, use u (0:n-2) */
    VSINT      n, /* dimensional of tridiagonal matrix */
    VScom   *com  /* common data */
)
{
    VSFLOAT a, b, c, delta, lb, tol1, tol2, break_tol,
         *gx, *gy, *gs, *e, *sc ;
    VSFLOAT p, q, r, s, scale, sigma, t, w, x, y, z0, z1 ;
    VSINT i, ip1, im1, k, l, l1, nm1, istart, nstart,
       *gi, *gj, *gk, *gf, *start_index ;
    short *gz ;
    VSINT ngivens,      /* number of fast Givens */
          nscale ;      /* number of scalings performed */
    int npass,          /* number of passes */
        max_npass ;
    VSparm *Parm = com->Parm ; /* parameters */
    VSgivens *G  = com->G ;    /* Givens rotations in diagonalization */

    /* output, not defined on input: */
    gx = G->gx ;      /* fast Givens x factor, VSFLOAT */
    gy = G->gy ;      /* fast Givens y factor, VSFLOAT */
    gz = G->gz ;      /* rotation type 1 or 2, short */
    gs = G->gs ;      /* scaling factors, VSFLOAT */
    gi = G->gi ;      /* row indices of scale factors, VSINT */
    gj = G->gj ;      /* indexed by QR pass, points into Givens factors, INT */
    gk = G->gk ;      /* index of last diagonal element in QR iteration, INT */
    gf = G->gf ;      /* index of first diagonal element in QR iteration */

    nm1 = n - 1 ;

    tol1 = Parm->qr_tol_rel ;  /* stop when u_{k-1}^2 <= tol1*|d_k| */
    tol2 = Parm->qr_tol_abs ;  /* stop when |u_{k-1}| <= tol2 assuming matrix
                                  scaled with absolute max element one */
    max_npass = VSMAX (30, Parm->qr_its*n) ;/* max QR algorithm passes */

    /* workspace vectors of size n for sc and start_index */
    sc = com->wx ;
    start_index = com->wi ;
    start_index [0] = 0 ;
    com->ngivens = ngivens = 0 ;
    com->npass = npass = 0 ;
    com->nscale = nscale = 0 ;

    /* ---------------------------------------------------------------------- */
    /* check for 1-by-1 case */
    /* ---------------------------------------------------------------------- */
    if (n == 1) return (0) ;

    /* ---------------------------------------------------------------------- */
    /* determine the problem scaling */
    /* ---------------------------------------------------------------------- */
    scale = vs_max(d, n) ;
    nstart = 0 ;
    for (i = 0 ; i < nm1 ; i++)
    {
        t = u [i] ;
        if ( scale < fabs (t) ) scale = fabs (t) ;
        if ( t == VSZERO )
        {
            nstart++ ;
            start_index [nstart] = i + 1 ;
        }
    }
    if (scale == 0 )
    {
        return (0) ;
    }
    scale = 1/scale ;
    /* ---------------------------------------------------------------------- */
    /* make a scaled copy of the input problem */
    /* ---------------------------------------------------------------------- */
    vs_scale (d, scale, n) ;
    vs_scale (u, scale, nm1) ;
    vs_initx (sc, VSONE, n) ;
    u [nm1] = 0 ; /* unused subdiagonal element */
    /* ---------------------------------------------------------------------- */
    /* inverse shifted QR method */
    /* ---------------------------------------------------------------------- */
    lb = Parm->qr_lower ; /* threshold in fast Givens for rescaling of rows */
    break_tol = lb*tol2 ;
    nscale = 1 ;
    gi [0] = 0 ;
    k = n ;
    l = nm1 ;

    while (k > 1)
    {
        l1 = 0 ;
        k = l ;
        l-- ;
        while ( k == start_index [nstart] )
        {
            k = l ;
            l-- ;
            nstart-- ;
            if ( k == 0 ) break ;
        }
        if ( k == 0 ) break ;

        while ( 1 )
        {
            if ( l1 >= 30 ) return (-5) ;
            t = sqrt(sc [l]*sc [k])*fabs(u [l]) ;
            if ( t <= tol2 ) break ;
            if ( t <= tol1*fabs(sc [k]*d [k]) ) break ;
            /* -------------------------------------------------------------- */
            /* log the start of a pass */
            /* -------------------------------------------------------------- */
            l1++ ;
            /* check if enough space was allocated for Givens rotations,
               if not, need to increase Parm->qr_its */
            if (npass >= max_npass) return (-6) ;
            gk [npass] = k ;
            istart = gf [npass] = start_index [nstart] ;
            gj [npass] = ngivens ;

            npass++ ;
            /* -------------------------------------------------------------- */
            /* compute shift sigma */
            /* -------------------------------------------------------------- */
            /* note that a, b, c, t, and delta */
            /* are not used outside this scope */
            a = d [k]*sc [k] ;
            b = d [l]*sc [l] ;
            c = sqrt(sc [l]*sc [k])*fabs(u [l]) ;
            delta = .5*(b-a) ;
            if (delta >= 0)
            {
                if (delta <= c)
                {
                    t = delta/c ;
                    sigma = a - c/(t+sqrt(t*t+1)) ;
                }
                else
                {
                    t = c/delta ;
                    sigma = a - c*t/(1+sqrt(t*t+1)) ;
                }
            }
            else
            {
                if (-delta <= c)
                {
                    t = -delta/c ;
                    sigma = a + c/(t+sqrt(t*t+1)) ;
                }
                else
                {
                    t = -c/delta ;
                    sigma = a + c*t/(1+sqrt(t*t+1)) ;
                }
            }

            /* -------------------------------------------------------------- */
            /* compute the initial p and q */
            /* -------------------------------------------------------------- */

            p = d [istart] - sigma/sc [istart] ;
            q = u [istart] ;
            /* -------------------------------------------------------------- */
            /* perform one pass */
            /* -------------------------------------------------------------- */
            for (i = istart; i < k; i++)
            { 
                /* test for break up of matrix */
                if ( i > istart )
                {
                    if ( fabs (u [i-1])*z0*z1 < break_tol )
                    {
                        nstart++ ;
                        start_index [nstart] = i ;
                    }
                }
	        /* ---------------------------------------------------------- */
	        /* scale the problem, if necessary */
	        /* ---------------------------------------------------------- */
                im1 = i - 1 ;
                z0 = sc [i] ;
                if ( z0 < lb )
                {
                    gi [nscale] = (i+1) ;
                    d [i] *= z0 ;
                    z0 = sqrt (z0) ;
                    gs [nscale] = z0 ;
                    nscale++ ;
                    p *= z0 ;
                    u [i] *= z0 ;
                    if ( i > 0 ) u [im1] *= z0 ;
                    z0 = VSONE ;
                    sc [i] = VSONE ;
                }
                ip1 = i + 1 ;
                z1 = sc [ip1] ;
                if ( z1 < lb )
                {
                    gi [nscale] = -(i+1) ;
                    d [ip1] *= z1 ;
                    z1 = sqrt (z1) ;
                    gs [nscale] = z1 ;
                    nscale++ ;
                    q *= z1 ;
                    u [i] *= z1 ;
                    u [ip1] *= z1 ;
                    if ( i > 0 ) a *= z1 ;
                    z1 = VSONE ;
                    sc [ip1] = VSONE ;
                }

                if (fabs (p) > fabs (q))
                {
                    r = z1/z0 ;
                    s = q/p ;
                    t = r*s*s ;

                    if (t <= 1)
                    {
                        t = 1/(1+t) ;
                        sc [i] = z0*t ;
                        sc [ip1]=z1*t ;
                        y = s ;
                        x = r*s ;
                        gx [ngivens] = x ;
                        gy [ngivens] = y ;
                        gz [ngivens] = 1 ;
                        ngivens++ ;

                        if (i > 0) u [im1] = u [im1] + x*a ;
                        t = d [i] ;
                        p = u [i] ;
                        q = t + x*p ;
                        r = p - y*t ;
                        w = d [ip1] - y*p ;
                        s = p + x*d [ip1] ;
                        u [i] = s - y*q ;
                        d [i] = q + x*s ;
                        d [ip1] = w - r*y ;
                        q = u [ip1] ;
                        a = x*q ;
                        p = w - sigma/z1 ;
                    }
                    else
                    {

                        t = t/(1+t) ;
                        r = z0/z1 ;
                        s = p/q ;
                        sc [i] = z1*t ;
                        sc [ip1] = z0*t ;
                        y = s ;
                        x = r*s ;
                        gx [ngivens] = x ;
                        gy [ngivens] = y ;
                        gz [ngivens] = 2 ;
                        ngivens++ ;

                        if (i > 0) u [im1] = x*u [im1] + a ;
                        t = d [i] ;
                        p = u [i] ;
                        q = x*t+p ;
                        r = y*p-t ;
                        w = y*d [ip1]-p ;
                        s = x*p+d [ip1] ;
                        u [i] = y*s-q ;
                        d [i] = x*q+s ;
                        d [ip1] = y*w-r ;
                        a = u [ip1] ;
                        q = a ;
                        u [ip1] = y*a ;
                        p = w - y*sigma/z1 ;
                    }

                }
                else
                {
                    if (q == 0)
                    {
                        gx [ngivens] = 0 ;
                        gy [ngivens] = 0 ;
                        gz [ngivens] = 1 ;
                        ngivens++ ;

                        p = d [ip1] - sigma/z1 ;
                        q = u [ip1] ;
                    }
                    else
                    {
                        r = z0/z1 ;
                        s = p/q ;
                        t = r*s*s ;

                        if (t < 1)
                        {
                            t = 1/(1+t) ;
                            sc [i] = z1*t ;
                            sc [ip1] = z0*t ;
                            y = s ;
                            x = r*s ;
                            gx [ngivens] = x ;
                            gy [ngivens] = y ;
                            gz [ngivens] = 2 ;
                            ngivens++ ;

                            if (i > 0) u [im1] = x*u [im1] + a ;
                            t = d [i] ;
                            p = u [i] ;
                            q = x*t+p ;
                            r = y*p-t ;
                            w = y*d [ip1] - p ;
                            s = x*p + d [ip1] ;
                            u [i] = y*s - q ;
                            d [i] = x*q + s ;
                            d [ip1] = y*w - r ;
                            a = u [ip1] ;
                            q = a ;
                            u [ip1] = y*a ;
                            p = w - y*sigma/z1 ;
                        }
                        else
                        {
                            t = t/(1+t) ;
                            r = z1/z0 ;
                            s = q/p ;
                            sc [i] = z0*t ;
                            sc [ip1] = z1*t ;
                            y = s ;
                            x = r*s ;
                            gx [ngivens] = x ;
                            gy [ngivens] = y ;
                            gz [ngivens] = 1 ;
                            ngivens++ ;

                            if (i > 0) u [im1] = u [im1] + x*a ;
                            t = d [i] ;
                            p = u [i] ;
                            q = t + x*p ;
                            r = p - y*t ;
                            w = d [ip1] - y*p ;
                            s = p + x*d [ip1] ;
                            u [i] = s - y*q ;
                            d [i] = q + x*s ;
                            d [ip1] = w - r*y ;
                            q = u [ip1] ;
                            a = x*q ;
                            p = w - sigma/z1 ;
                        }
                    }
                }
            }
            gi [nscale] = 0 ;
            nscale++ ;
        }
    }

    /* remove the initial scaling applied to the tridiagonal matrix */
    scale = 1./scale ;
    for (i = 0 ; i < n ; i++ )
    {
        d [i] *= sc [i]*scale ;
    }
    com->nscale = nscale ;

    if ( ngivens > 0 ) 
    {
        for (i = 0 ; i < n ; i++)
        {
            gs [nscale+i] = sqrt (sc [i]) ;
        }
    }

    com->ngivens = ngivens ;
    com->npass = npass ;
    return (0) ;
}

/* ==========================================================================
   === vs_DenseMult =========================================================
   ==========================================================================

   Compute y = Vx where V is m by n */

void vs_DenseMult
(
    VSFLOAT *y,     /* m by 1 product, output */
    VSFLOAT *x,     /* n by 1 given vector */
    VSFLOAT *V,     /* dense m by n matrix */
    VSINT    m,     /* number of rows */
    VSINT    n      /* number of columns */
)
{
    VSINT i, j ;
    VSFLOAT *Vp, t ;

    if ( n < 1 ) return ;
    t = x [0] ;
    for (i = 0; i < m; i++) y [i] = V [i]*t ;
    Vp = V+m ;
    for (j = 1; j < n; j++)
    {
        t = x [j] ;
        for (i = 0; i < m; i++) y [i] += Vp [i]*t ;
        Vp += m ;
    }
}

/* ==========================================================================
   === vs_tDenseMult ==========================================================
   ==========================================================================

   Compute y' = x'V where V is m by n */

void vs_tDenseMult
(
    VSFLOAT *y,     /* n by 1 product, output */
    VSFLOAT *x,     /* m by 1 given vector */
    VSFLOAT *V,     /* dense m by n matrix */
    VSINT    m,     /* number of rows */
    VSINT    n      /* number of columns */
)
{
    VSINT i, j ;
    VSFLOAT *Vp, t ;

    Vp = V ;
    for (j = 0; j < n; j++)
    {
        t = 0. ;
        for (i = 0; i < m; i++)
        {
            t += Vp [i]*x [i] ;
        }
        y [j] = t ;
        Vp += m ;
    }
}

/* ========================================================================== */
/* === vs_GivensMult ======================================================== */
/* ========================================================================== */

/* Multiply row vector x on the left by the Givens rotations generated
   during the QR algorithm.  The rotations are applied in the same order
   that they were generated during the QR algorithm. */
    
void vs_GivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT       n  /* dimension of matrix */
)
{   
    VSINT i, j, k, f, *gk, *gf, *gi, ngivens, nscale, pass ;
    short *gz, *Gz ;
    VSFLOAT xxt, xxs, *gx, *gy, *gs, *Gx, *Gy, *Gs ;

    /* ---------------------------------------------------------------------- */
    /* Read in the needed arrays */
    /* ---------------------------------------------------------------------- */

    /* input */
    if ( npass == 0 ) return ; /* diagonal matrix */
    gx = G->gx ;       /* fast Givens x factor */
    gy = G->gy ;       /* fast Givens y factor */
    gz = G->gz ;       /* rotation type 1 or 2 */
    gs = G->gs ;       /* scaling factors */
    gi = G->gi ;       /* row indices for scaling factors */
    gk = G->gk ;       /* index of last diagonal element in QR iteration */
    gf = G->gf ;       /* index of first diagonal element in QR iteration */

    ngivens = 0 ;
    nscale = 0 ;
    for (pass = 0; pass < npass; pass++)
    {
        k = gk [pass] ;
        f = gf [pass] ;

        /* apply Givens rotation */
        Gx = gx+(ngivens-f) ;
        Gy = gy+(ngivens-f) ;
        Gz = gz+(ngivens-f) ;
        nscale++ ;
        i = abs (gi [nscale]) - 1 ;
        for (j = gf [pass]; j < k; j++)
        {
            if ( i == j )
            {
                if ( gi [nscale] > 0 )
                {
                    x [j] *= gs [nscale] ;
                    nscale++ ;
                    i = abs (gi [nscale]) - 1 ;
                    if ( i == j )
                    {
                        x [j+1] *= gs [nscale] ;
                        nscale++ ;
                        i = abs (gi [nscale]) - 1 ;
                    }
                }
                else
                {
                    x [j+1] *= gs [nscale] ;
                    nscale++ ;
                    i = abs (gi [nscale]) - 1 ;
                }
            }
            xxt = x [j] ;
            xxs = x [j+1] ;
            if ( Gz [j] == 1 )
            {
                x [j]   = xxt + Gx [j] * xxs ;
                x [j+1] = xxs - Gy [j] * xxt ;
            }
            else
            {
                x [j]   = Gx [j] * xxt + xxs ;
                x [j+1] = Gy [j] * xxs - xxt ;
            }
        }
        ngivens += (k-f) ;
    }

    /* final scaling */
    nscale++ ;
    Gs = gs+nscale ;
    for (i = 0; i < n ; i++) x [i] *= Gs [i] ;
}

/* ==========================================================================
   === vs_tGivensMult =======================================================
   ==========================================================================
   Multiply a column vector x on the right by the Givens rotations generated
   during the QR algorithm.  The rotations are applied in the reverse order
   that they were generated during the QR algorithm. */

void vs_tGivensMult
(
    VSFLOAT    *x, /* the vector to which the rotations are applied */
    VSgivens   *G, /* Givens structure for storing Given rotation info */
    VSINT   npass, /* number of QR iterations used to compute eigenvalues */
    VSINT  nscale, /* number of scaling operations performed */
    VSINT       n  /* dimension of matrix */
)
{
    VSINT jj, *gj ;
    VSINT i, j, k, kk, kf, *gi, *gk, *gf ;
    short *gz, *Gz ;
    VSFLOAT *gx, *gy, *gs, *Gx, *Gy, *Gs ;

    /* input */
    if ( npass == 0 ) return ; /* diagonal matrix */
    gx = G->gx ;       /* fast Givens x factor */
    gy = G->gy ;       /* fast Givens y factor */
    gz = G->gz ;       /* rotation type 1 or 2 */
    gs = G->gs ;       /* scaling factors */
    gi = G->gi ;       /* indexed by QR pass, points into scale factors */
    gj = G->gj ;       /* indexed by QR pass, points into Givens factors */
    gk = G->gk ;       /* index of last diagonal element in QR iteration */
    gf = G->gf ;       /* index of first diagonal element in QR iteration */

    Gs = gs+nscale ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] *= Gs [j] ;
    }
    nscale-- ;

    for (k = npass-1 ; k >= 0 ; k--)
    {
        jj = gj [k] ;
        kk = gk [k] -1 ;
        kf = gf [k] ;
        Gx = gx+(jj-kf) ;
        Gy = gy+(jj-kf) ;
        Gz = gz+(jj-kf) ;
        nscale-- ;
        i = abs (gi [nscale]) - 1 ;

        for (j = kk; j >= kf; j--)
        {
            double xt = x [j] ;
            double xs = x [j+1] ;

            if (Gz [j] == 1)
            {
                x [j]   = xt - Gy [j] * xs ;
    	        x [j+1] = xs + Gx [j] * xt ;
            }
            else
            {
                x [j]   = Gx [j] * xt - xs ;
                x [j+1] = Gy [j] * xs + xt ;
            }
            /* scale */
            if ( j == i )
            {
                if ( gi [nscale] < 0 )
                {
                    x [j+1] *= gs [nscale] ;
                    nscale-- ;
                    i = abs (gi [nscale]) - 1 ;
                    if ( i == j )
                    {
                        x [j] *= gs [nscale] ;
                        nscale-- ;
                        i = abs (gi [nscale]) - 1 ;
                    }
                }
                else
                {
                    x [j] *= gs [nscale] ;
                    nscale-- ;
                    i = abs (gi [nscale]) - 1 ;
                }
            }
        }
    }
}

/* ==========================================================================
   === vs_printA ============================================================
   ==========================================================================
    Print a sparse matrix
   ========================================================================== */
void vs_printA
(
    VSINT   const     n, /* number of cols in A */
    VSLINT  const   *Ap, /* size n+1, column pointers */
    VSINT   const   *Ai, /* size Ap [n], row indices for A */
    VSFLOAT const   *Ax, /* size Ap [n], numerical entries of A */
    int     const  diag, /* 1 = add diag */
    char    const *what  /* name of the matrix */
)
{
    VSINT j ;
    VSLINT p, q, t ;
    VSFLOAT *colmax ;
    colmax = malloc ((n+1)*sizeof (VSFLOAT)) ;
    printf ("Temp = [\n") ;
    p = 0 ;
    for (j = 0; j < n; j++)
    {
        q = Ap [j+1] ;
        if ( Ax != NULL )
        {
            t = 0 ;
            for (; p < q; p++)
            {
                printf("%ld %ld %e\n",
                        (long) (Ai [p]+1), (long) (j+1), Ax [p]) ;
                if ( Ax [p] > t ) t = Ax [p] ;
            }
            colmax [j+1] = t ;
        }
        else
        {
            for (; p < q; p++)
            {
                printf ("%ld %ld %ld\n",
                        (long) (Ai [p]+1), (long) (j+1), (long) -1) ;
            }
            colmax [j+1] = 1 ;
        }
    }
    if ( diag )
    {
        for (j = 1; j <= n; j++) printf ("%i %i %e\n", j, j, colmax [j]) ;
    }
    printf ("] ;\n") ;
    printf ("%s = sparse (Temp(:, 1), Temp(:, 2), Temp(:, 3)) ;\n", what) ;
}

void vs_printxMATLAB
(
    VSFLOAT *x, /* numerical entries in the vector */
    VSINT    n, /* dimension of the vector */
    char *what  /* name of the vector */
)
{
    VSINT j ;
    printf ("%s = [\n", what) ;
    for (j = 0; j < n; j++)
    {
        printf ("%26.16e\n", x [j]) ;
    }
    printf ("] ;\n") ;
}

/* napheap codes */
/*
   =========================================================================
   ============================== NAPHEAP ==================================
   =========================================================================

       ________________________________________________________________
      |Solve a separable convex quadratic knapsack problem of the form:|
      |                                                                |
      | min .5x'Dx - c'x  subject to  lo <= x <= hi, blo <= a'x <= bhi |
      |                                                                |
      |       where lo and hi are vectors, blo and bhi are scalars,    |
      |       and D is a diagonal matrix with nonnegative diagonal.    |
      |                                                                |
      |       Copyright November 1, 2019 by                            |
      |   Timothy A. Davis, William W. Hager, and James T. Hungerford  |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/Software        |
      |                                                                |
      |  Disclaimer: The views expressed are those of the authors and  |
      |              do not reflect the official policy or position of |
      |              the Department of Defense or the U.S. Government. |
      |                                                                |
      |      Approved for Public Release, Distribution Unlimited       |
      |________________________________________________________________|

    The constraints are rewritten as

               lo <= x <= hi,    blo <= b <= bhi,    a'x = b

    The code solves the dual problem obtained by introducing
    a multiplier lambda for the constraint a'x = b.  The dual function is

    L(lambda) = inf {.5x'Dx - c'x + lambda (a'x - b):
                                 lo <= x <= hi, blo <= b <= bhi}

    The dual function is concave, and the solution to the original (primal)
    problem can be constructed from the maximizer of the dual function.
    If blo < bhi, then there is a discontinuity at lambda = 0.  D is a
    diagonal matrix, and if d = diag(D) contains zeros, there could be
    additional points of discontinuity.  The algorithm iterates on lambda
    until 0 lies in the subdifferential of L(lambda).  The optimal x is
    constructed from the minimizers of the dual function.  Three algorithms
    are implemented: a break point searching algorithm, a variable fixing
    algorithm, and a Newton/secant algorithm.  A starting guess for the
    optimal multiplier can be provided.  If the user does not have a good
    starting guess, at most Parm->K (default 20) iterations of a Newton-type
    algorithm are used to generate a starting guess. The Newton-type method
    terminates when the root is bracketed.

    By default, Newton's method is used for the startup (Parm->newton =
    TRUE).  Set Parm->newton to FALSE to choose the variable fixing
    algorithm for startup. After the startup iterations, the breakpoint
    algorithm is used to compute the final solution. When Newton's method is
    used for the startup, all entries of d must be positive (otherwise an error
    is returned), so to solve problems with zeros in d, you must set
    Parm->newton to FALSE.  No value of d can be negative, for any method.

    The user can either provide work arrays or let the code malloc work arrays.
    The necessary size for the work arrays depends on the choice of the
    algorithm, as can be seen in the memory allocation routine.  The most
    memory that the code will require is an NAPINT work array of size 4*n+2
    and a real work array of size 5*n where n is the problem dimension.

    The default parameter values used by the software are given in the
    napheap_default routine.  To alter the default parameter values, the
    user can first initialize the parameters to their default values using
    napheap_default, make any changes to the values in the parameter
    structure, and provide the updated parameter structure as an input to the
    napheap code. Note that TRUE = 1 and FALSE = 0.

    The code records the number of iterations for either the variable fixing or
    Newton method, and the number of break points that are passed through
    during the solution process.  These statistics are returned when the user
    provides a nonnull structure for the statistics argument of the code.  The
    statistics and parameters can be printed with napheap_print_stat and
    napheap_print_parms.

    If the user wishes to override the default definition of infinity
    given in Include/napheap.h, modify the OPTFLAGS as explained in
    Lib/napheap.mk.
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|

      Alternative licenses are available. Contact William Hager for details.

      Reference:
         T. A. Davis, W. W. Hager, and J. T. Hungerford, An Efficient
         Hybrid Algorithm for the Separable Convex Quadratic Knapsack Problem,
         see Doc/napheap_paper.pdf
*/

/* ========================================================================== */
/* === user-callable functions ============================================== */
/* ========================================================================== */

/* There are only four user-callable functions:
    napheap_default  - sets default parameters for napheap
    napheap          - solves a separable convex quadratic knapsack problem
    napheap_print    - print napheap parameters and statistics
    naptest          - exhaustive test of napheap code (in Test directory)
*/

/* ========================================================================== */
/* === napheap  ============================================================= */
/* ========================================================================== */

int napheap  /* Return status of solution process. To determine the status, use
                napheap_print_status (napdata) */
(
    NAPdata *napdata
)
{
    int    status ;
    VSINT it, k, nk, nuf, nkf ;

    status = NAPHEAP_STATUS_OK ;   /* default status */
    NAPINT const n = napdata->n ;  /* problem dimension */
    if ( n <= 0 )
    {
        status = NAPHEAP_STATUS_INVALID_N ;
        napdata->Parm->PrintStat = FALSE ;
        napdata->Stat->kerror = n ;
        return (nap_wrapup (status, napdata, napdata->lambda)) ;
    }

    /* if y is NULL, then it is treated as zero */
    if ( napdata->c == NULL )
    {
        napdata->c = (NAPFLOAT *) napheap_malloc (&status, n, sizeof(NAPFLOAT));
        if ( status )
        {
            status = NAPHEAP_STATUS_OUT_OF_MEMORY ;
            napdata->Parm->PrintStat = FALSE ;
            return (nap_wrapup (status, napdata, napdata->lambda)) ;
        }
        napdata->c_created = napdata->c ;
        for (k = 0; k < n; k++) napdata->c [k] = NAPZERO ;
    }

    /* create the solution array if the user did not provide one */
    if ( napdata->x == NULL )
    {
        napdata->x = (NAPFLOAT *) napheap_malloc (&status, n, sizeof(NAPFLOAT));
        if ( status )
        {
            status = NAPHEAP_STATUS_OUT_OF_MEMORY ;
            napdata->Parm->PrintStat = FALSE ;
            return (nap_wrapup (status, napdata, napdata->lambda)) ;
        }
        napdata->x_created = napdata->x ;
    }

    /* grab the problem data */
    NAPFLOAT        *x = napdata->x ;      /* solution */
    NAPFLOAT    lambda = napdata->lambda ; /* linear constraint multiplier */
    /* In the code, y is used for linear cost vector. For consistency with the
       rest of SuiteOPT, c is used for the linear cost vector in the input data
       structure. */
    NAPFLOAT const  *y = napdata->c ;      /* linear cost vector */
    NAPFLOAT const  *d = napdata->d ;      /* diagonal of objective Hessian */
    NAPFLOAT const  *a = napdata->a ;      /* linear constraint vector */
    NAPFLOAT const *lo = napdata->lo ;     /* lower bounds for x */
    NAPFLOAT const *hi = napdata->hi ;     /* upper bounds for x */
    NAPstat      *Stat = napdata->Stat ;   /* structure for statistics */
    NAPparm      *Parm = napdata->Parm ;   /* structure for parameters */
    NAPFLOAT    *xWork = napdata->xWork ;  /* NULL = code creates NAPFLOAT mem*/
    NAPINT      *iWork = napdata->iWork ;  /* NULL = code creates NAPINT mem */

    int const a_is_one = (a == NULL) ? TRUE : FALSE ;
    /* NOTE: napdata->blo and napdata->bhi are extracted later,
             after napheap_check */

    /* data from a prior run if it exists */
    NAPFLOAT const akakdk = napdata->akakdk ;

    /* internal variable used in the code */
    int use_newton ; /* TRUE if starting guess is Newton iteration */
    int goright ;    /* TRUE if optimal multiplier to right of current lambda */
    int crossed ;    /* TRUE when iterate crosses root */
    /* Rout and Lout are only used when a diagonal element vanishes */
    int Rout ;       /* becomes FALSE when optimal lambda < lambdaR and hence,
                        lambdaR cannot be optimal */
    int Lout ;       /* becomes FALSE when optimal lambda > lambdaL and hence,
                        lambdaL cannot be optimal */
    int K ;          /* number of iterations in Newton's or variable fixing
                        methods */
    int k_gt_1 ;     /* TRUE if K > 1 */
    int recompute_slope ;     /* TRUE if slope at current lambda not known */
    NAPFLOAT maxax, minax ;
    NAPFLOAT *ada = NULL ;      /* size n, ada[j] == aj*aj */
    NAPFLOAT *ady = NULL ;      /* size n, ady[j] == aj*yj */
    NAPFLOAT *ad = NULL ;       /* size n, ad[j] = aj/dj */
    NAPFLOAT *ay = NULL ;       /* size n, ay[j] = yj/aj */
    /* breakpts [j] is defined only when j is in the list of undetermined
       variables.  It is the first break point between the current lambda
       and the root */
    NAPFLOAT *breakpts = NULL ; /* size n */
    NAPFLOAT *breaknext = NULL ;/* size n, stores next break point */
    NAPFLOAT *br_hi = NULL ;    /* size n. br_hi[j] hi break pt for x[j] */
    NAPFLOAT *br_lo = NULL ;    /* size n, br_lo[j] lo break pt for x[j] */
    NAPFLOAT b ;                /* a'x = b */
    NAPFLOAT domlo ;            /* lower bound for domain of dual function */
    NAPFLOAT domhi ;            /* upper bound for domain of dual function */
    NAPFLOAT slope ;            /* slope at lambda, exclude subdiff. term */
    NAPFLOAT slopelo ;          /* lower subdifferential at lambda */
    NAPFLOAT slopehi ;          /* upper subdifferential at lambda */
    NAPFLOAT slope0 ;           /* slope at lambda = 0 excluding subdiff term */
    NAPFLOAT slope0lo ;         /* lower subdifferential at lambda = 0 */
    NAPFLOAT slope0hi ;         /* upper subdifferential at lambda = 0 */
    NAPFLOAT prevlambda ;       /* value of lambda in previous iteration */
    NAPFLOAT pprevlambda ;      /* value of lambda in pre-previous iteration */
    NAPFLOAT lambdaL ;          /* lower bound for optimal multiplier */
    NAPFLOAT lambdaR ;          /* upper bound for optimal multiplier */
    NAPFLOAT dL ;               /* slope at lambdaL */
    NAPFLOAT dR ;               /* slope at lambdaR */
    NAPFLOAT prevslope ;        /* slope in the previous iteration */
    NAPFLOAT pprevslope ;       /* slope in the pre-previous iteration */
    int  prevstep;              /* type of step at each iteration */
    NAPFLOAT num ;              /* numerator of variable fixing ratio */
    NAPFLOAT den ;              /* denominator of variable fixing ratio */
    NAPFLOAT numFV ;            /* initial num of variable fixing ratio */
    NAPFLOAT denFV ;            /* initial den of variable fixing ratio */

    NAPINT *uf = NULL ;         /* size n, unfixed variables */
    NAPINT *kf = NULL ;         /* size n, known free variables */
    NAPINT *free_heap = NULL ;  /* size n+1, heap of free variables */
    NAPINT *bound_heap = NULL ; /* size n+1, heap of bound variables */
    NAPINT *iwork ;             /* internal NAPINT workspace */
    NAPFLOAT *xw2 = NULL ;      /* internal NAPFLOAT workspace */

    int start_napsearch;        /* TRUE when break point search used */
    int start_napsearch_next ;  /* TRUE if napsearch starts next iteration */

    NAPINT nnewton ;            /* number of newton steps taken */
    NAPINT nsecant ;            /* number of secant steps taken */
    NAPINT nvarfix ;            /* number of vf steps taken during loop
                                   iterations (these are slow steps) */
    NAPINT nbound ;             /* number of indices in bound heap */
    NAPINT nfree ;              /* number of indices in free heap */
    /* use_lambda is TRUE if the user provides a starting guess for the
       multiplier associated with the constraint blo <= a'x <= bhi */
//printf ("lambda input: %e\n", lambda) ;
    int use_lambda = ( lambda >= NAPINF ) ? FALSE : TRUE ;
//printf("use_lambda: %i\n", use_lambda) ;

    /* ---------------------------------------------------------------------- */
    /*  Initialization and checks */
    /* ---------------------------------------------------------------------- */
    if ( Parm == NULL )
    {
        status = NAPHEAP_STATUS_MISSING_PARM ;
        return (nap_wrapup (status, napdata, lambda)) ;
    }
    if ( Stat == NULL )
    {
        status = NAPHEAP_STATUS_MISSING_STAT ;
        return (nap_wrapup (status, napdata, lambda)) ;
    }
    if ( Parm->PrintParm )
    {
        napheap_print_parm (napdata) ;
    }

    int const use_prior_data = Parm->use_prior_data ;

    /* when using prior data, add statistics to prior totals,
       when not using prior data, initialize statistics to 0 */
    if ( !use_prior_data )
    {
        Stat->kerror  = 0 ;
        Stat->lobad   = 0 ;
        Stat->hibad   = 0 ;
        Stat->dbad    = 0 ;
        Stat->nkf     = 0 ;
        Stat->nfree   = 0 ;
        Stat->nbound  = 0 ;
        Stat->nbrks   = 0 ;
        Stat->nrefine = 0 ;
        Stat->nvarfix = 0 ;
        Stat->nnewton = 0 ;
        Stat->nsecant = 0 ;
    }

    /* d_is_zero is TRUE if entries in d are all zero (input d is ignored).
       if d is NULL, then it is treated as identically zero. */
    const int d_is_zero = ((d == NULL) && (Parm->d_is_one == FALSE)) ?
                          TRUE : Parm->d_is_zero ;

    /* d_is_one is TRUE if entries in d are all one (input d is ignored) */
    const int d_is_one = Parm->d_is_one ;

    /* check for consistency in the parameter values */
    if ( (d_is_zero == TRUE) && (d_is_one == TRUE) )
    {
        status = NAPHEAP_STATUS_CONFLICTING_PARM_D ;
        return (nap_wrapup (status, napdata, lambda)) ;
    }

    /* If lo == NULL, then loExists is always FALSE (lo = -INF), otherwise
       loExists = Parm->loExists */
    const int loExists = (lo == NULL) ? FALSE : Parm->loExists ;

    /* If hi == NULL, then hiExists is always FALSE (hi = INF), otherwise
       hiExists = Parm->hiExists */
    const int hiExists = (hi == NULL) ? FALSE : Parm->hiExists ;

    const int lohiExist = TRUE ;

    /* d_is_pos TRUE if d is identically 1 or all entries in d are positive.
       When d_is_pos is FALSE, there could be both zero and positive entries
       in d. */
    if ( d_is_zero == TRUE ) k = FALSE ;
    else                     k = Parm->d_is_pos ;
    const int d_is_pos = (d_is_one) ? TRUE : k ;

    /* check bounds and d if Parm->check is TRUE */
    if ( Parm->check )
    {
        /* Note: if Parm->check is false, then blo is permitted to be larger
           than bhi,  and we silently swap the two bounds.  However, when
           Parm->check is true, this case generates an error, since it might
           be a mistake on the part of the user.  */
        if ( napdata->blo > napdata->bhi )
        {
            status = NAPHEAP_STATUS_INVALID_B ;
            return (nap_wrapup (status, napdata, lambda)) ;
        }
        for (k = 0; k < n; k++)
        {
            if ( d_is_pos || !d_is_zero )
            {
                if ( (d [k] < NAPZERO) || ((d [k] == NAPZERO) && d_is_pos) )
                {
                    status = NAPHEAP_STATUS_INVALID_D ;
                    Stat->kerror = k ;
                    Stat->dbad = d [k] ;
                    return (nap_wrapup (status, napdata, lambda)) ;
                }
            }
        }
    }

    /* switch bounds if napdata->blo > napdata->bhi */
    NAPFLOAT const blo =
                   (napdata->blo > napdata->bhi) ? napdata->bhi : napdata->blo ;
    NAPFLOAT const bhi =
                   (napdata->blo > napdata->bhi) ? napdata->blo : napdata->bhi ;

    status = NAPHEAP_STATUS_OK ;   /* default status */

    /* ---------------------------------------------------------------------- */
    /* check for basic errors */
    /* ---------------------------------------------------------------------- */
    if (n == 0)
    {
        /* empty problem, do nothing */
        return (nap_wrapup (status, napdata, lambda)) ;
    }

    /* Newton's method requires all (d>0). !d_is_pos is an error when
       Newton's method is specified. */
    const int newton = (!d_is_pos) ? FALSE : Parm->newton ;

    K = Parm->K ;
    if ( K < 0 ) K = 0 ; /* silently treat this as K = 0 */

    /* ---------------------------------------------------------------------- */
    /* initialize parameters */
    /* ---------------------------------------------------------------------- */

    /* use_lambda is TRUE if starting lambda is provided */
    if ( use_lambda && (lambda == NAPZERO) && (blo < bhi) )
    {
        use_lambda = FALSE ; /* slope at lambda = 0 is always evaluated */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate NAPINT workspace */
    /* ---------------------------------------------------------------------- */

    if ( use_prior_data == TRUE )
    {
        iwork = iWork ;
    }
    else if ( iWork == NULL ) /* create memory if it is not provided by user */
    {
        iwork = (NAPINT *) napheap_malloc (&status, (4*n+2), sizeof (NAPINT)) ;
        /* Need space for uf, kf, free_heap, and bound_heap.  uf has size n
           while the heaps have size n+1: heap [0], heap [1], ... , heap [n] */
        if ( status )
        {
            return (nap_wrapup (NAPHEAP_STATUS_OUT_OF_MEMORY, napdata, lambda));
        }
        napdata->iWork_created = napdata->iWork = iwork ;
    }
    else
    {
        /* use the user-provided workspace, iWork */
        iwork = iWork ;
    }

    /* assign pointers to NAPINT work arrays */
    uf = iwork ; iwork += n ;               /* uf, size n */
    kf = iwork ; iwork += n ;               /* kf, size n */
    bound_heap = iwork ; iwork += (n+1) ;   /* bound_heap, size n+1 */
    free_heap = iwork ;  iwork += (n+1) ;   /* free_heap, size n+1 */


    /* ---------------------------------------------------------------------- */
    /* allocate NAPFLOAT workspace */
    /* ---------------------------------------------------------------------- */

    if ( (K == 0) || d_is_pos )
    {
        if ( d_is_one )
        {
            /* 3n workspace */
            status = nap_work (n, xWork, &xw2,
                NULL, &ay, NULL,  NULL, &breakpts, &breaknext) ;
        }
        else
        {
            /* 4n workspace */
            status = nap_work (n, xWork, &xw2,
                &ad, &ay, NULL,  NULL, &breakpts, &breaknext) ;
        }
    }
    else if ( !d_is_zero )
    {
        /* 5n workspace */
        status = nap_work (n, xWork, &xw2,
            &ada, &ady, &br_hi, &br_lo, &breakpts, NULL);
    }
    else if ( K > 1 ) /* d_is_zero and K > 1 */
    {
        /* 3n workspace */
        status = nap_work (n, xWork, &xw2,
            &ada, &ady, NULL,  NULL, &breakpts, NULL) ;
    }
    else    /* d_is_zero and K = 1 */
    {
        /* 1n workspace */
        status = nap_work (n, xWork, &xw2,
            NULL, NULL, NULL,  NULL, &breakpts, NULL) ;
    }
    if (status != NAPHEAP_STATUS_OK)
    {
        return (nap_wrapup (status, napdata, lambda)) ;
    }
    /* if nap_work created some memory, then save pointer in napdata */
    if ( xw2 != NULL )
    {
        napdata->xWork_created = napdata->xWork = xw2 ;
    }

    /* ---------------------------------------------------------------------- */
    /* remaining initializations */
    /* ---------------------------------------------------------------------- */

    /* initialize variables for feasibility testing */
    maxax = NAPZERO ;   /* stores hi range of possible values for a'x */
    minax = NAPZERO ;   /* stores lo range of possible values for a'x */

    domhi = NAPINF ;    /* upper bound on domain where dual max achieved */
    domlo = -NAPINF ;   /* lower bound on domain where dual max achieved */

    /* initialize other variables */
    slope = NAPZERO ;   /* slope of dual func ignoring break points at lambda */
    slope0 = NAPZERO ;  /* slope of the dual function at lambda = 0 */
    b = NAPZERO ;       /* adjustments to right side due to fixed variables */
    lambdaL = +NAPINF ; /* most negative break point */
    lambdaR = -NAPINF ; /* most positive break point */
    dL = +NAPINF ;      /* value of slope at lambdaL */
    dR = -NAPINF ;      /* value of slope at lambdaR */
    num = NAPZERO ;     /* numerator of lambda estimate in variable fixing */
    den = NAPZERO ;     /* denominator of lambda estimate in variable fixing */
    nuf = 0 ;           /* number of unfixed variables */
    nkf = 0 ;           /* number of known free variables */
    nbound = 0 ;        /* no variables known to be at bounds */
    nfree = 0 ;         /* no variables known to be free */

    /* ====================================================================== */
    /* === solve the knapsack problem (3 cases: d > 0, d = 0, complement) === */
    /* ====================================================================== */

    /* lo and/or hi could have infinite valued components. */

    if ( d_is_pos )
    {
        /* ================================================================== */
        /* === case when (d > 0) ============================================ */
        /* ================================================================== */



        /* if no starting guess given, compute one by variable fixing */
        if ( !use_lambda )
        {

//printf("no start guess given, use variable fixing, do not use lambda\n");

            if ( (use_prior_data == TRUE) && (akakdk < NAPINF) )
            {
                den = akakdk ;
                if ( a_is_one )
                {
                    num = 0 ;
                    for (k = 0; k < n; k++)
                    {
                        num += y [k] ;
                    }
                }
                else
                {
                    num = nap_dot (y, a, n) ;
                }
//printf("den: %e num: %e b: %e\n", den, num, b) ;
                if ( blo == bhi )
                {
                    b = blo ;
                }
                else /* blo < bhi */
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( yk < VSZERO ) continue ;
                        if ( a_is_one )
                        {
                            if ( yk > VSONE )
                            {
                                slope0++ ;
                            }
                            else
                            {
                                slope0 += yk ;
                            }
                        }
                        else
                        {
                            if ( yk > VSONE )
                            {
                                slope0 += a [k] ;
                            }
                            else
                            {
                                slope0 += a [k]*yk ;
                            }
                        }
                    }
                    /* determine bracketting interval and b */
                    BRACKET ;
                }
            }

            else /* starting from scratch */
            {
                if ( blo == bhi )
                {
                    for (k = 0; k < n; k++)
                    {
                        if ( a_is_one )
                        {
                            den++ ;
                            num +=  y [k] ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            den += ak*ak ;
                            num +=  y [k]*ak ;
                        }
                    }
                    b = blo ;
                    if ( Parm->return_data )
                    {
                        napdata->akakdk = den ;
                    }
                }
                else /* blo < bhi */
                {
//printf ("start from scratch\n") ;
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            den++ ;
                            num += yk ;
                            if ( yk < VSZERO )
                            {
                                continue ;
                            }
                            if ( yk > VSONE )
                            {
                                slope0++ ;
                            }
                            else
                            {
                                slope0 += yk ;
                            }
                        }
                        else /* ak = a [k] */
                        {
                            const NAPFLOAT ak = a [k] ;
                            den += ak*ak ;
                            num +=  yk*ak ;
                            if ( yk < VSZERO )
                            {
                                continue ;
                            }
                            if ( yk > VSONE )
                            {
                                slope0 += ak ;
                            }
                            else
                            {
                                slope0 += ak*yk ;
                            }
                        }
                    }
                    if ( Parm->return_data )
                    {
                        napdata->akakdk = den ;
                    }
//printf ("do not use lambda\n") ;
//printf ("slope0: %e slope: %e den: %e\n", slope0, slope, den) ;
                    /* determine bracketting interval and b */
                    BRACKET ;
                }
             
                if ( Parm->return_data )
                {
                    napdata->akakdk = den ;
                }
            }
            lambda = (num - b)/den ;
//printf ("b: %e num: %e lambda: %e\n", b, num, lambda) ;
            numFV = num ; /* numerator for variable fixing */
            denFV = den ; /* denominator for variable fixing */

            /* compute slope at starting point, slope at lambda = 0,
               and den for Newton */
            den = NAPZERO ;
            if ( blo == bhi )
            {
                for (k = 0; k < n; k++)
                {
                    const NAPFLOAT yk = y [k] ;
                    if ( a_is_one )
                    {
                        const NAPFLOAT t = yk - lambda ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope++ ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += t ;
                            den++ ;
                        }
                        if ( yk > VSONE )
                        {
                            slope0++ ;
                        }
                        else if ( yk > VSZERO )
                        {
                            slope0 += yk ;
                        }
                    }
                    else
                    {
                        const NAPFLOAT ak = a [k] ;
                        const NAPFLOAT t = yk - lambda*ak ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope += ak ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += ak*t ;
                            den += ak*ak ;
                        }
                        if ( yk > VSONE )
                        {
                            slope0 += ak ;
                        }
                        else if ( yk > VSZERO )
                        {
                            slope0 += ak*yk ;
                        }
                    }
                }
//printf ("slope: %e slope0: %e\n", slope, slope0) ;
                /* determine bracketting interval */
                BRACKET ;
            }
            else /* blo < bhi */
            {
                for (k = 0; k < n; k++)
                {
                    if ( a_is_one )
                    {
                        const NAPFLOAT t = y [k] - lambda ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope++ ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += t ;
                            den++ ;
                        }
                    }
                    else
                    {
                        const NAPFLOAT ak = a [k] ;
                        const NAPFLOAT t = y [k] - ak*lambda ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope += ak ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += ak*t ;
                            den += ak*ak ;
                        }
                    }
                }
            }
        }
        else /* use lambda, compute slope at starting guess and at lambda = 0*/
        {
//printf("use lambda: %e\n", lambda) ;
            /* no need to compute slope at lambda = 0 when blo = bhi */
            if ( blo == bhi )
            {
                for (k = 0; k < n; k++)
                {
                    if ( a_is_one )
                    {
                        const NAPFLOAT t = y [k] - lambda ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                            continue ;
                        }
                        if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope++ ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += t ;
                            den++ ;
                        }
                    }
                    else
                    {
                        const NAPFLOAT ak = a [k] ;
                        const NAPFLOAT t = y [k] - lambda*ak ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                            continue ;
                        }
                        if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope += ak ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += ak*t ;
                            den += ak*ak ;
                        }
                    }
                }
                b = blo ;
                lambdaL =-NAPINF ;
                lambdaR = NAPINF ;
                dL = NAPINF ;
                dR =-NAPINF ;
            }
            else /* blo < bhi, also compute slope at lambda = 0 */
            {
                for (k = 0; k < n; k++)
                {
                    const NAPFLOAT yk = y [k] ;
                    if ( a_is_one )
                    {
                        const NAPFLOAT t = yk - lambda ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else
                        {
                            if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope++ ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += t ;
                                den++ ;
                            }
                        }
                        if ( yk < VSZERO )
                        {
                            continue ;
                        }
                        if ( yk > VSONE )
                        {
                            slope0++ ;
                        }
                        else
                        {
                            slope0 += yk ;
                        }
                    }
                    else /* ak != 1 */
                    {
                        const NAPFLOAT ak = a [k] ;
                        const NAPFLOAT t = yk - lambda*ak ;
                        if ( t < VSZERO )
                        {
                            x [k] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [k] = VSONE ;
                            slope += ak ;
                        }
                        else
                        {
                            x [k] = t ;
                            slope += ak*t ;
                            den += ak*ak ;
                        }
                        if ( yk < VSZERO )
                        {
                            continue ;
                        }
                        if ( yk > VSONE )
                        {
                            slope0 += ak ;
                        }
                        else
                        {
                            slope0 += ak*yk ;
                        }
                    }
                }
                /* determine bracketting interval and b */
                BRACKET ;
            }
//printf ("slope: %25.15e slope0: %25.15e den: %25.15e\n", slope, slope0, den) ;
        }
//printf ("b: %25.15e\n", b) ;
        slope -= b ;
//printf ("new slope: %25.15e\n", slope) ;

        /* update bracketting interval based on slope at starting guess */
        if ( slope > NAPZERO )
        {
            if ( dL < slope )
            {
                lambda = lambdaL ;
                slope = dL ;
            }
            else
            {
                lambdaL = lambda ;
                dL = slope ;
            }
        }
        else if ( slope < NAPZERO )
        {
            if ( dR > slope )
            {
                lambda = lambdaR ;
                slope = dR ;
            }
            else
            {
                lambdaR = lambda ;
                dR = slope ;
            }
        }
        else /* given lambda is optimal */
        {
            status = NAPHEAP_STATUS_OK ;
            return (nap_wrapup (status, napdata, lambda)) ;
        }
//printf ("slope0: %e\n", slope0) ;
//printf ("slope: %e\n", slope) ;
//printf ("lambda: %e\n", lambda) ;
//printf ("K = %i\n", K) ;

        /* At this point, we have a starting guess that was either given
           or was computed by a variable fixing iteration. x has been
           set in accordance with the bounds at the starting guess.
           The slope is known at both the starting guess and at lambda = 0
           (unless blo = bhi and use_lambda = TRUE), and we have an
           interval that brackets the root. We now fix bound variables,
           remove components of the a array that vanish, initialize the
           uf array of unfixed indices, and compute either a Newton or
           a variable fixing iteration. */

        if ( K == 0 )
        {
            /* x values are set according to the starting guess. If
               lambda = 0 is the current iterate, then recompute x using
               lambda = 0 */
            if ( lambda == NAPZERO )
            {
                for (k = 0; k < n; k++)
                {
                    const NAPFLOAT t = y [k] ;
                    if ( t < VSZERO )
                    {
                        x [k] = VSZERO ;
                    }
                    else if ( t > VSONE )
                    {
                        x [k] = VSONE ;
                    }
                    else
                    {
                        x [k] = t ;
                    }
                    ay [k] = (a_is_one) ? y [k] : y [k]/a [k] ;
                    uf [nuf] = k ;
                    nuf++ ;
                }
            }
            else /* lambda != 0 */
            {
                for (k = 0; k < n; k++)
                {
                    ay [k] = (a_is_one) ? y [k] : y [k]/a [k] ;
                    uf [nuf] = k ;
                    nuf++ ;
                }
            }
            status = napsearch (x, &lambda, lambdaL, lambdaR, b, slope,
                                ay, ad, a, lo, hi, n, breakpts, breaknext,
                                nuf, &nkf, &nfree, &nbound, uf, kf,
                                free_heap, bound_heap, Stat, Parm,
                                loExists, hiExists) ;

            if ( status == NAPHEAP_STATUS_OK )
            {
                napsolution1 (x, y, d, a, lo, hi, lambda, slope, kf, nkf,
                              free_heap, bound_heap, nfree, nbound, d_is_one) ;
            }

            return (nap_wrapup (status, napdata, lambda)) ;
        }

        use_newton = FALSE ;
        prevlambda = lambda ;
        prevslope = slope ;
/*printf ("newton: %i\n", newton) ;*/
        if ( newton )
        {
            if ( den != NAPZERO ) /* note: lambda=0 => using den of prevlambda*/
            {
                lambda += Parm->newton_scale*slope/den ; /* scaled Newton */
/*printf ("new lambda: %e lambdaL: %e lambdaR: %e\n",
           lambda, lambdaL, lambdaR) ; */
                if ( (lambda >= lambdaL) && (lambda <= lambdaR) )
                {
                    use_newton = TRUE ;
                }
            }

            /* try secant if newton failed */
            if ( !use_newton )
            {
                if ( (lambdaL > -NAPINF) && (lambdaR < NAPINF) ) /* use secant*/
                {
                    /* express the secant update relative to the side with
                       slope of smallest magnitude */
                    if ( dL < -dR )
                    {
                        lambda = lambdaL - dL*(lambdaL - lambdaR)/(dL-dR) ;
                    }
                    else
                    {
                        lambda = lambdaR - dR*(lambdaL - lambdaR)/(dL-dR) ;
                    }
                    use_newton = TRUE ;
                }
            }
        }

        int zero_bracket = FALSE ;
        if ( (lambdaL == NAPZERO) || (lambdaR == NAPZERO) )
        {
            zero_bracket = TRUE ; /* one side of bracketing interval is 0 */
        }
/* printf ("zero_bracket: %i\n", zero_bracket) ;
printf ("fix variables and evaluate next iterate\n") ;*/

        /* fix variables and evaluate next iterate */
        slope = NAPZERO ;
        if ( use_newton )
        {
            den = NAPZERO ;
            if ( !zero_bracket ) /* zero not in bracketting interval */
            {
                if ( prevslope > NAPZERO )
                {
                    for (k = 0; k < n; k++)
                    {
                        if ( a_is_one )
                        {
                            if ( x [k] == VSZERO )
                            {
                                continue ;
                            }
                            const NAPFLOAT yk = y [k] ;
                            ay [k] = yk ;
                            const NAPFLOAT t = yk - lambda ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else
                            {
                                if ( t > VSONE )
                                {
                                    x [k] = VSONE ;
                                    slope += VSONE ;
                                }
                                else
                                {
                                    x [k] = t ;
                                    slope += t ;
                                    den += VSONE ;
                                }
                            }
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSZERO )
                            {
                                continue ;
                            }
                            const NAPFLOAT yk = y [k] ;
                            ay [k] = yk/ak ;
                            const NAPFLOAT t = yk - lambda*ak ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else
                            {
                                if ( t > VSONE )
                                {
                                    x [k] = VSONE ;
                                    slope += ak ;
                                }
                                else
                                {
                                    x [k] = t ;
                                    slope += ak*t ;
                                    den += ak*ak ;
                                }
                            }
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
                else /* prevslope < 0 */
                {
                    for (k = 0; k < n; k++)
                    {
                        if ( a_is_one )
                        {
                            if ( x [k] == VSONE )
                            {
                                b -= VSONE ;
                                continue ;
                            }
                            const NAPFLOAT yk = y [k] ;
                            ay [k] = yk ;
                            const NAPFLOAT t = yk - lambda ;
                            if ( t  > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += VSONE ;
                            }
                            else if ( t < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += t ;
                                den += VSONE ;
                            }
                        }
                        else
                        {
                            if ( x [k] == VSONE )
                            {
                                b -= a [k] ;
                                continue ;
                            }
                            const NAPFLOAT ak = a [k] ;
                            const NAPFLOAT yk = y [k] ;
                            ay [k] = yk/ak ;
                            const NAPFLOAT t = yk - lambda*ak ;
                            if ( t  > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += ak ;
                            }
                            else if ( t < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += ak*t ;
                                den += ak*ak ;
                            }
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
            }
            else if ( prevlambda == NAPZERO ) /*bracket is [0, INF)|(-INF, 0]*/
            {
                /* fix variables at lambda = 0 */
                if ( prevslope > NAPZERO )
                {
                    for (k = 0; k < n; k++)
                    {
                        if ( a_is_one )
                        {
                            const NAPFLOAT yk = y [k] ;
                            if ( yk <= VSZERO )
                            {
                                x [k] = VSZERO ;
                                continue ;
                            }
                            ay [k] = yk ;
                            const NAPFLOAT t = yk - lambda ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += VSONE ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += t ;
                                den += VSONE ;
                            }
                        }
                        else
                        {
                            const NAPFLOAT yk = y [k] ;
                            if ( yk <= VSZERO )
                            {
                                x [k] = VSZERO ;
                                continue ;
                            }
                            const NAPFLOAT ak = a [k] ;
                            ay [k] = yk/ak ;
                            const NAPFLOAT t = yk - lambda*ak ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += ak ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += ak*t ;
                                den += ak*ak ;
                            }
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
                else /* prevslope < 0 */
                {
                        for (k = 0; k < n; k++)
                        {
                            const NAPFLOAT yk = y [k] ;
                            if ( a_is_one )
                            {
                                if ( yk >= VSONE )
                                {
                                    b -= VSONE ;
                                    x [k] = VSONE ;
                                    continue ;
                                }
                                ay [k] = yk ;
                                const NAPFLOAT t = yk - lambda ;
                                if ( t  > VSONE )
                                {
                                    x [k] = VSONE ;
                                    slope += VSONE ;
                                }
                                else if ( t < VSZERO )
                                {
                                    x [k] = VSZERO ;
                                }
                                else
                                {
                                    x [k] = t ;
                                    slope += t ;
                                    den += VSONE ;
                                }
                            }
                            else
                            {
                                const NAPFLOAT ak = a [k] ;
                                if ( yk >= VSONE )
                                {
                                    b -= ak ;
                                    x [k] = VSONE ;
                                    continue ;
                                }
                                ay [k] = yk/ak ;
                                const NAPFLOAT t = yk - lambda*ak ;
                                if ( t  > VSONE )
                                {
                                    x [k] = VSONE ;
                                    slope += ak ;
                                }
                                else if ( t < VSZERO )
                                {
                                    x [k] = VSZERO ;
                                }
                                else
                                {
                                    x [k] = t ;
                                    slope += ak*t ;
                                    den += ak*ak ;
                                }
                            }
                            uf [nuf] = k ;
                            nuf++ ;
                        }
                }
            }
            else /* bracket is [0, lambda] or [lambda, 0] */
            {
//printf ("bracket is [0, lambda] or [lambda, 0]\n") ;
                if ( prevslope > NAPZERO )
                {
//printf ("prevslope > 0\n") ;
                    for (k = 0; k < n; k++)
                    {
                        if ( x [k] == VSZERO )
                        {
                            continue ;
                        }
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( yk > VSONE )
                            {
                                x [k] = VSONE ;
                                b -= VSONE ;
                                continue ;
                            }
                            ay [k] = yk ;
                            const NAPFLOAT t = yk - lambda ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += VSONE ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += t ;
                                den += VSONE ;
                            }
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( yk > VSONE )
                            {
                                x [k] = VSONE ;
                                b -= ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                            const NAPFLOAT t = yk - lambda*ak ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += ak ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += ak*t ;
                                den += ak*ak ;
                            }
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
//printf("slope: %e den: %e\n", slope, den) ;
                }
                else /* prevslope < 0 */
                {
                    for (k = 0; k < n; k++)
                    {
                        if ( a_is_one )
                        {
                            if ( x [k] == VSONE )
                            {
                                b -= VSONE ;
                                continue ;
                            }
                            const NAPFLOAT yk = y [k] ;
                            if ( yk < VSZERO )
                            {
                                x [k] = VSZERO ;
                                continue ;
                            }
                            ay [k] = yk ;
                            const NAPFLOAT t = yk - lambda ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += VSONE ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += t ;
                                den += VSONE ;
                            }
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSONE )
                            {
                                b -= ak ;
                                continue ;
                            }
                            const NAPFLOAT yk = y [k] ;
                            if ( yk < VSZERO )
                            {
                                x [k] = VSZERO ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                            const NAPFLOAT t = yk - lambda*ak ;
                            if ( t  < VSZERO )
                            {
                                x [k] = VSZERO ;
                            }
                            else if ( t > VSONE )
                            {
                                x [k] = VSONE ;
                                slope += ak ;
                            }
                            else
                            {
                                x [k] = t ;
                                slope += ak*t ;
                                den += ak*ak;
                            }
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
            }
            /* If nuf = 0, then all the variables are fixed.
               Optimal solution has been found. */
            if ( nuf == 0 )
            {
                status = NAPHEAP_STATUS_OK ;
                lambda = NAPZERO ;
                return (nap_wrapup (status, napdata, lambda)) ;
            }
            slope -= b ;
        }
        else /* variable fixing iterate */
        {
            /* need to compute num and dem if use_lambda is TRUE */
            if ( use_lambda )
            {
                den = NAPZERO ;
                num = NAPZERO ;
                for (k = 0; k < n; k++)
                {
                    if ( a_is_one )
                    {
                        den += VSONE ;
                        num += y [k] ;
                    }
                    else
                    {
                        const NAPFLOAT ak = a [k] ;
                        den += ak*ak ;
                        num += y [k]*ak ;
                    }
                }
            }
            else /* num and den already computed above */
            {
                num = numFV ;
                den = denFV ;
            }
            /* this differs from Newton both in the formula for num and den
               and due to the fact that there is no lambda in the bracketting
               interval */
            if ( !zero_bracket ) /* zero not in bracketting interval */
            {
                if ( prevslope > NAPZERO )
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( x [k] == VSZERO )
                            {
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            ay [k] = yk ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSZERO )
                            {
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
                else /* prevslope < 0 */
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( x [k] == VSONE )
                            {
                                b -= VSONE ;
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            ay [k] = yk ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSONE )
                            {
                                b -= ak ;
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
            }
            else if ( prevlambda == NAPZERO ) /*bracket is [0, INF) |(-INF, 0]*/
            {
                /* fix variables at lambda = 0 */
                if ( prevslope > NAPZERO )
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( yk <= VSZERO )
                            {
                                x [k] = VSZERO ;
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            ay [k] = yk ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( yk <= VSZERO )
                            {
                                x [k] = VSZERO ;
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
                else /* prevslope < 0 */
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( yk >= VSONE )
                            {
                                b -= VSONE ;
                                x [k] = VSONE ;
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            ay [k] = yk ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( yk >= VSONE )
                            {
                                b -= ak ;
                                x [k] = VSONE ;
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
            }
            else /* bracket is [0, lambda] or [lambda, 0] */
            {
                if ( prevslope > NAPZERO )
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( x [k] == VSZERO )
                            {
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            if ( yk > VSONE )
                            {
                                x [k] = VSONE ;
                                b -= VSONE ;
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            ay [k] = yk ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSZERO )
                            {
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            if ( yk > VSONE )
                            {
                                x [k] = VSONE ;
                                b -= ak*VSONE ;
                                den -= ak*ak ;
                                num -= yk*ak ;
                                continue ;
                            }
                            ay [k] = yk/ak ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
                else /* prevslope < 0 */
                {
                    for (k = 0; k < n; k++)
                    {
                        const NAPFLOAT yk = y [k] ;
                        if ( a_is_one )
                        {
                            if ( x [k] == VSONE )
                            {
                                b -= VSONE ;
                                den -= VSONE ;
                                num -= yk ;
                                continue ;
                            }
                            const NAPFLOAT r = y [k] ;
                            if ( r < VSZERO )
                            {
                                x [k] = VSZERO ;
                                den -= VSONE ;
                                num -= r ;
                                continue ;
                            }
                            ay [k] = r ;
                        }
                        else
                        {
                            const NAPFLOAT ak = a [k] ;
                            if ( x [k] == VSONE )
                            {
                                b -= ak ;
                                den -= ak*ak ;
                                num -= y [k]*ak ;
                                continue ;
                            }
                            const NAPFLOAT r = y [k]/ak ;
                            if ( r < VSZERO )
                            {
                                x [k] = VSZERO ;
                                den -= ak*ak ;
                                num -= y [k]*ak ;
                                continue ;
                            }
                            ay [k] = r ;
                        }
                        uf [nuf] = k ;
                        nuf++ ;
                    }
                }
            }
            /* If nuf = 0, then all the variables are fixed. Either the
               problem is infeasible or the optimal solution has been found. */
            if ( nuf == 0 )
            {
                status = NAPHEAP_STATUS_OK ;
                lambda = NAPZERO ;
                return (nap_wrapup (status, napdata, lambda)) ;
            }
            if ( newton )
            {
                lambda = (num - b)/den ;
                den = NAPZERO ;
                LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                if ( a_is_one )
                {
                    const NAPFLOAT t = y [j] - lambda ;
                    if ( t < VSZERO )
                    {
                        x [j] = VSZERO ;
                    }
                    else if ( t > VSONE )
                    {
                        x [j] = VSONE ;
                        slope += VSONE ;
                    }
                    else
                    {
                        x [j] = t ;
                        slope += t ;
                        den += VSONE ;
                    }
                }
                else
                {
                    const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
                    const NAPFLOAT t = y [j] - lambda*aj ;
                    if ( t < VSZERO )
                    {
                        x [j] = VSZERO ;
                    }
                    else if ( t > VSONE )
                    {
                        x [j] = VSONE ;
                        slope += aj ;
                    }
                    else
                    {
                        x [j] = t ;
                        slope += aj*t ;
                        den += aj*aj ;
                    }
                }
                END_LOOP ;
                slope -= b ;
            }
        }

        /* ================================================================== */
        /* Begin newton algorithm */
        /* ================================================================== */
        if ( newton )
        {
            /* update bracketting interval based on slope at current iterate */
            if ( slope > NAPZERO )
            {
                lambdaL = lambda ;
                dL = slope ;
            }
            else if ( slope < NAPZERO )
            {
                lambdaR = lambda ;
                dR = slope ;
            }
            else /* given lambda is optimal */
            {
                status = NAPHEAP_STATUS_OK ;
                return (nap_wrapup (status, napdata, lambda)) ;
            }

            /* We now have an iterate lambda, the slope and den at lambda,
               the array uf of unfixed indices */

            /* ============================================================== */
            /* Start Newton method */
            /* ============================================================== */
            /* Perform up to K iterations of Newton's method. If we cannot
               perform a Newton step since L'' vanishes, then perform a
               secant step as long as both sides of the bracketing interval
               are finite. If one side is infinite, then perform a variable
               fixing iteration. Switch to the breakpoint algorithm if
               no variables are fixed in an iteration, the iteration limit K
               is reached, or either a Newton or a secant iteration is
               performed which generates a point on the opposite side of
               the root. If a variable fixing step takes us to the other
               side of the root, then we try to perform one additional Newton
               iteration before switching to the breakpoint algorithm.
               When the Newton iteration fails to fix a variable,
               then we switch to the break point algorithm. */
            status = NAPHEAP_STATUS_KEEP_LOOKING ; /* default status */

            /* number of steps */
            nnewton    = 0 ;
            nsecant    = 0 ;
            nvarfix    = 0 ;
            prevstep = NONE ;

            /* When lambda crosses root, start_napsearch becomes true. */
            start_napsearch = FALSE ;

            /* If a variable fixing iteration crosses the root, then
               start_napsearch_next becomes TRUE */
            start_napsearch_next = FALSE ;

            for ( it = 1; it <= K; it++ )
            {
                crossed = FALSE ; /* changes to TRUE when root is crossed */
                if ( prevslope*slope < NAPZERO )
                {
                    crossed = TRUE ;
                }

                /* ========================================================== */
                /* Check to see if we should switch to napsearch */
                /* ========================================================== */
                if ( start_napsearch_next || (it == K) )
                {
                    start_napsearch = TRUE;
                }
                else if ( (it > 1) ||
                          (use_lambda && (fabs (slope) <= fabs (prevslope))) )
                {
                    if ( crossed && (prevstep != VARFIX) )
                    {
                        start_napsearch = TRUE;
                    }
                    else if ( crossed ) /* start napsearch next iteration */
                    {
                        start_napsearch_next = TRUE;
                    }
                }
                pprevslope = prevslope ;
                prevslope = slope;    /* save current slope */
                pprevlambda = prevlambda ;
                prevlambda = lambda ; /* save current lambda */

                /* ========================================================== */
                /* Step computation */
                /* ========================================================== */
                /* If L'' = 0, use secant if both sides of the bracketing
                   interval are finite otherwise use variable fixing. */
                if ( den == NAPZERO )
                {
                    /* Take a variable fixing step if bracket interval inf.
                       This will generate a point inside [lambdaL, lambdaR]
                       since we also fix variables. */
                    if ( lambdaL == -NAPINF || lambdaR == NAPINF )
                    {
                        prevstep = VARFIX ;
                        nvarfix++ ;
                        den = NAPZERO ;
                        num = NAPZERO ;
                        if ( d_is_one )
                        {
                            if ( slope > NAPZERO )
                            {
                                LOOP_OVER_UNFIXED ;
                                {
                                    if ( a_is_one )
                                    {
                                        if ( x [j] == VSZERO )
                                        {
                                            continue ;
                                        }
                                        den += VSONE ;
                                        num += y [j] ;
                                    }
                                    else
                                    {
                                        const NAPFLOAT aj = (a_is_one) ?
                                                                 VSONE : a [j] ;
                                        if ( x [j] == VSZERO )
                                        {
                                            continue ;
                                        }
                                        den += aj*aj ;
                                        num += aj*y [j] ;
                                    }
                                }
                                END_UNFIXED_LOOP ;
                            }
                            else /* slope < 0 */
                            {
                                LOOP_OVER_UNFIXED ;
                                {
                                    if ( a_is_one )
                                    {
                                        if ( x [j] == VSONE )
                                        {
                                            b-- ;
                                            continue;
                                        }
                                        den++ ;
                                        num += y [j] ;
                                    }
                                    else
                                    {
                                        const NAPFLOAT aj = a [j] ;
                                        if ( x [j] == VSONE )
                                        {
                                            b -= aj ;
                                            continue;
                                        }
                                        den += aj*aj ;
                                        num += aj*y [j] ;
                                    }
                                }
                                END_UNFIXED_LOOP ;
                            }
                        }
                        else
                        {
                            if ( slope > NAPZERO )
                            {
                                LOOP_OVER_UNFIXED ;
                                {
                                    if ( a_is_one )
                                    {
                                        if ( x [j] == VSZERO )
                                        {
                                            continue ;
                                        }
                                        den++ ;
                                        num += ay [j] ;
                                    }
                                    else
                                    {
                                        const NAPFLOAT aj = a [j] ;
                                        if ( x [j] == VSZERO )
                                        {
                                            continue ;
                                        }
                                        const NAPFLOAT t = aj*aj ;
                                        den += t ;
                                        num += t*ay [j] ;
                                    }
                                }
                                END_UNFIXED_LOOP ;
                            }
                            else /* slope < 0 */
                            {
                                LOOP_OVER_UNFIXED ;
                                {
                                    if ( a_is_one )
                                    {
                                        if ( x [j] == VSONE )
                                        {
                                            b -= VSONE ;
                                            continue;
                                        }
                                        den++ ;
                                        num += ay [j] ;
                                    }
                                    else
                                    {
                                        const NAPFLOAT aj = a [j] ;
                                        if ( x [j] == VSONE )
                                        {
                                            b -= aj ;
                                            continue;
                                        }
                                        const NAPFLOAT t = aj*aj ;
                                        den += t ;
                                        num += t*ay [j] ;
                                    }
                                }
                                END_UNFIXED_LOOP ;
                            }
                        }

                        if ( nuf == 0 )
                        {
                            /* all components fixed, but slope != 0 */
                            status = NAPHEAP_STATUS_INFEASIBLE ;
                            return (nap_wrapup (status, napdata, lambda)) ;
                        }
                        /* den != 0 since k > -1 */
                        lambda = (num - b)/den ;
                    }
                    else /* den = 0 and finite interval, perform secant step */
                    {
                        if ( prevstep == SECANT1 )
                        {
                            prevstep = SECANT2 ;
                        }
                        else
                        {
                            prevstep = SECANT1 ;
                        }
                    }
                }
                else /*den != 0*/
                {
                    lambda += Parm->newton_scale*slope/den; /* scaled Newton */
                    if ( (lambda <= lambdaR) && (lambda >= lambdaL) )
                    {
                        prevstep = NEWTON ;
                        nnewton++;
                    }
                    else
                    {
                        if ( prevstep == SECANT1 )
                        {
                            prevstep = SECANT2 ;
                        }
                        else
                        {
                            prevstep = SECANT1 ;
                        }
                    }

                }

                /* Use secant step if lambda lands outside (lambdaL, lambdaR)
                   or den = 0.  This can only happen when lambda is
                   generated by a Newton iteration and the interval is finite.
                   We use two different secant steps. The first is based on
                   the function values and slopes at the ends of the bracketting
                   interval. The second is based on the last two iterates.
                   We alternate between these two ways to compute the secant
                   step. */
                if ( prevstep == SECANT2 )
                {
                    nsecant++ ;

                    /* express the secant update relative to the side with
                       slope of smallest magnitude */
                    const NAPFLOAT t = prevslope - pprevslope ;
                    if ( t != NAPZERO )
                    {
                        if ( fabs (prevslope) < fabs (pprevslope) )
                        {
                            lambda = prevlambda - prevslope*
                                 (prevlambda - pprevlambda)/t ;
                        }
                        else
                        {
                            lambda = pprevlambda - pprevslope*
                                 (prevlambda - pprevlambda)/t ;
                        }
                    }
                    /* if the step was no good, then switch to a secant step
                       based on the bracketting interval */
                    if ((t == NAPZERO)||(lambda < lambdaL)||(lambda > lambdaR))
                    {
                        nsecant-- ;
                        prevstep = SECANT1 ;
                    }
                }
                if ( prevstep == SECANT1 )
                {
                    nsecant++ ;

                    /* express the secant update relative to the side with
                       slope of smallest magnitude */
                    if ( dL < -dR )
                    {
                        lambda = lambdaL - dL*(lambdaL - lambdaR)/(dL-dR) ;
                    }
                    else
                    {
                        lambda = lambdaR - dR*(lambdaL - lambdaR)/(dL-dR) ;
                    }
                }
                /* lambda has now been determined (either by Newton, secant,
                   or variable fixing) */

                /* ========================================================== */
                /* Fix variables using previous lambda and slope              */
                /* evaluate L' and L'' at the new lambda                      */
                /* ========================================================== */

                slope = NAPZERO ;
                den   = NAPZERO ;
                /* If the previous step was a variable fixing iteration,
                   then we skip the fixing of variables since this was done
                   while the iterate was computed. */
                if ( prevstep == VARFIX )
                {
                    if ( start_napsearch )
                    {
                        if ( d_is_one )
                        {
                            /* fixing done already, compute slope, set x */
                            LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                                if ( a_is_one )
                                {
                                    const NAPFLOAT t = y [j] - lambda ;
                                    if ( t < VSZERO )
                                    {
                                        x [j] = VSZERO ;
                                    }
                                    else
                                    {
                                        if ( t > VSONE )
                                        {
                                            x [j] = VSONE ;
                                            slope += VSONE ;
                                        }
                                        else
                                        {
                                            x [j] = t ;
                                            slope += t ;
                                        }
                                    }
                                }
                                else
                                {
                                    const NAPFLOAT aj = (a_is_one) ?
                                                                 VSONE : a [j] ;
                                    const NAPFLOAT t = y [j] - lambda*aj ;
                                    if ( t < VSZERO )
                                    {
                                        x [j] = VSZERO ;
                                    }
                                    else
                                    {
                                        if ( t > VSONE )
                                        {
                                            x [j] = VSONE ;
                                            slope += aj ;
                                        }
                                        else
                                        {
                                            x [j] = t ;
                                            slope += aj*t ;
                                        }
                                    }
                                }
                            END_LOOP ;
                        }
                        else
                        {
                            /* fixing done already, compute slope, set x */
                            LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                                if ( a_is_one )
                                {
                                    const NAPFLOAT t = ay [j] - lambda ;
                                    if ( t < VSZERO )
                                    {
                                        x [j] = VSZERO ;
                                    }
                                    else
                                    {
                                        if ( t > VSONE )
                                        {
                                            x [j] = VSONE ;
                                            slope += VSONE ;
                                        }
                                        else
                                        {
                                            x [j] = t ;
                                            slope += t ;
                                        }
                                    }
                                }
                                else
                                {
                                    const NAPFLOAT aj = a [j] ;
                                    const NAPFLOAT t = (ay [j] - lambda)*aj ;
                                    if ( t < VSZERO )
                                    {
                                        x [j] = VSZERO ;
                                    }
                                    else
                                    {
                                        if ( t > VSONE )
                                        {
                                            x [j] = VSONE ;
                                            slope += aj ;
                                        }
                                        else
                                        {
                                            x [j] = t ;
                                            slope += aj*t ;
                                        }
                                    }
                                }
                            END_LOOP ;
                        }
                    }
                    else  /* start_napsearch = FALSE, prevstep = VARFIX, compute
                             slope & den for Newton step, fixing done earlier */
                    {
                        LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                            if ( a_is_one )
                            {
                                const NAPFLOAT t = y [j] - lambda ;
                                if ( t < VSZERO )
                                {
                                    x [j] = VSZERO ;
                                }
                                else
                                {
                                    if ( t > VSONE )
                                    {
                                        x [j] = VSONE ;
                                        slope += VSONE ;
                                    }
                                    else
                                    {
                                        x [j] = t ;
                                        slope += t ;
                                        den += VSONE ;
                                    }
                                }
                            }
                            else
                            {
                                const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
                                const NAPFLOAT t = y [j] - lambda*aj ;
                                if ( t < VSZERO )
                                {
                                    x [j] = VSZERO ;
                                }
                                else
                                {
                                    if ( t > VSONE )
                                    {
                                        x [j] = VSONE ;
                                        slope += aj ;
                                    }
                                    else
                                    {
                                        x [j] = t ;
                                        slope += aj*t ;
                                        den += aj*aj ;
                                    }
                                }
                            }
                        END_LOOP ;
                    }
                }
                else if ( prevslope > NAPZERO )
                { /* previous step was Newton, optimal multiplier to right */
                    if ( start_napsearch )
                    { /* fix variables, no need to compute den */
                        LOOP_OVER_UNFIXED ;
                        {
                            const NAPFLOAT aj = ( a_is_one ) ? VSONE : a [j] ;
                            FIX_LO ;
                            SLOPE_UPDATE_LO ;
                        }
                        END_UNFIXED_LOOP ;
                    }
                    else /* start_napsearch = FALSE, slope > 0 */
                    { /* continue Newton iteration, fix variables, compute den*/
                        LOOP_OVER_UNFIXED ;
                        {
                            const NAPFLOAT aj = ( a_is_one ) ? VSONE : a [j] ;
                            FIX_LO ;
                            SLOPE_UPDATE_LO_DEN ;
                        }
                        END_UNFIXED_LOOP ;
                    }
                }
                else /* prevslope < 0 */
                { /* Newton step previously, optimal multiplier to the left */
                    if ( start_napsearch )
                    { /* fix variables and compute break points */
                        LOOP_OVER_UNFIXED ;
                        {
                            const NAPFLOAT aj = ( a_is_one ) ? VSONE : a [j] ;
                            FIX_HI ;
                            SLOPE_UPDATE_HI ;
                        }
                        END_UNFIXED_LOOP ;
                    }
                    else /* start_napsearch = FALSE, slope < 0 */
                    { /* continue Newton iteration, fix variables, compute den*/
                        LOOP_OVER_UNFIXED ;
                        {
                            const NAPFLOAT aj = ( a_is_one ) ? VSONE : a [j] ;
                            FIX_HI ;
                            /* update slope and compute den*/
                            SLOPE_UPDATE_HI_DEN ;
                        }
                        END_UNFIXED_LOOP ;
                    }
                }
                slope -= b ;

                /* update bracketting interval */
                if ( slope > NAPZERO ) /* optimal mult. to the right */
                {
                    lambdaL = lambda ;
                    dL = slope ;
                }
                else if ( slope < NAPZERO ) /* optimal mult. is to the left */
                {
                    lambdaR = lambda ;
                    dR = slope ;
                }
                else /* slope = 0 and lambda is optimal */
                {
                    status = NAPHEAP_STATUS_OK ;
                    /* compute optimal solution for unfixed variables */
                    LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                        const NAPFLOAT aj = ( a_is_one ) ? VSONE : a [j] ;
                        NAPFLOAT t = y [j] - lambda*aj ;
                        if ( t < VSZERO )
                        {
                            t = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            t = VSONE ;
                        }
                        x [j] = t ;
                    END_LOOP ;
                    break ;
                }

                /* all variables fixed and slope != 0 */
                if ( nuf == 0 )
                {
                    status = NAPHEAP_STATUS_INFEASIBLE ;
                    return (nap_wrapup (status, napdata, lambda)) ;
                }

                /* When no variables are fixed, stop the Newton iteration
                   and start napsearch. */
                if ( start_napsearch || (nuf == nk) )
                {
                    break ;
                }
            }
            Stat->nnewton = nnewton ;
            Stat->nsecant = nsecant ;
            Stat->nvarfix = nvarfix;
        }
        else /* variable fixing algorithm */
        {
            prevslope = NAPZERO ;  /* dummy value, do not exit 1st iteration*/

            /* Perform up to K iterations of the variable fixing algorithm */
            for (it = 1; it <= K; it++)
            {
                /* evaluate new lambda and the slope */
                lambda = (num-b)/den ;
                crossed = FALSE ;
                if ( slope*prevslope < NAPZERO )
                {
                    crossed = TRUE ;
                }
                prevslope = slope ;
                slope = NAPZERO ;
                LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
                    if ( a_is_one )
                    {
                        const NAPFLOAT t = y [j] - lambda ;
                        if ( t < VSZERO )
                        {
                            x [j] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [j] = VSONE ;
                            slope += VSONE ;
                        }
                        else
                        {
                            slope += t ;
                            x [j] = t ;
                        }
                    }
                    else
                    {
                        const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
                        const NAPFLOAT t = y [j] - lambda*aj ;
                        if ( t < VSZERO )
                        {
                            x [j] = VSZERO ;
                        }
                        else if ( t > VSONE )
                        {
                            x [j] = VSONE ;
                            slope += aj ;
                        }
                        else
                        {
                            slope += aj*t ;
                            x [j] = t ;
                        }
                    }
                END_LOOP ;
                slope -= b ;

                if ( slope > NAPZERO )     /* optimal mult. is to the right */
                {
                    if ( lambdaL < lambda )
                    {
                        lambdaL = lambda ;
                        dL = slope ;
                    }
                    else
                    {
                        lambda = lambdaL ;
                        slope = dL ;
                    }

                    if ( it == K )
                    {
                        if ( nuf == 0 ) /* optimal solution found */
                        {
                            status = NAPHEAP_STATUS_OK ;
                            break ;
                        }
                        status = NAPHEAP_STATUS_KEEP_LOOKING ;
                        break ;
                    }
                    /* fix variables and compute next step */
                    LOOP_OVER_UNFIXED ;
                    {
                        const NAPFLOAT aj = (a_is_one) ? VSONE : a [j];
                        if ( aj > NAPZERO )
                        {
                            VAR_FIX_LO ;
                        }
                        else /* aj < NAPZERO */
                        {
                            VAR_FIX_HI ;
                        }
                    }
                    END_UNFIXED_LOOP ;
                }
                else if ( slope < NAPZERO ) /* optimal mult. is to the left */
                {
                    if ( lambdaR > lambda )
                    {
                        lambdaR = lambda ;
                        dR = slope ;
                    }
                    else
                    {
                        lambda = lambdaR ;
                        slope = dR ;
                    }
                    if ( it == K )
                    {
                        status = NAPHEAP_STATUS_KEEP_LOOKING ;
                        break ;
                    }

                    /* fix variables and compute next step */
                    LOOP_OVER_UNFIXED ;
                    {
                        const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
                        if ( aj > NAPZERO )
                        {
                            VAR_FIX_HI ;
                        }
                        else /* aj < NAPZERO */
                        {
                            VAR_FIX_LO ;
                        }
                    }
                    END_UNFIXED_LOOP ;
                }
                else /* optimal solution found (slope is zero) */
                {
                    status = NAPHEAP_STATUS_OK ;
                    break ;
                }

                /* all variables fixed and slope != 0 */
                if ( (nuf == 0) && (nkf == 0) )
                {
                    status = NAPHEAP_STATUS_INFEASIBLE ;
                    return (nap_wrapup (status, napdata, lambda)) ;
                }

                if ( nuf == nk ) /* no variables were fixed */
                {
                    /* optimal solution found */
                    status = NAPHEAP_STATUS_OK ;
                    break ;
                }
            }
            Stat->nvarfix = it ;
        }

        /* if optimal solution not found, perform breakpoint search */
        if ( status == NAPHEAP_STATUS_KEEP_LOOKING )
        {
            status = napsearch (x, &lambda, lambdaL, lambdaR, b, slope,
                                ay, ad, a, lo, hi, n, breakpts, breaknext,
                                nuf, &nkf, &nfree, &nbound, uf, kf,
                                free_heap, bound_heap, Stat, Parm,
                                loExists, hiExists) ;

            if ( status == NAPHEAP_STATUS_OK )
            {
                napsolution1 (x, y, d, a, lo, hi, lambda, slope, kf, nkf,
                     free_heap, bound_heap, nfree, nbound, d_is_one) ;
            }
        }
        else /* optimal solution found */
        {
            /* heaps not built, set any unfixed variables and free variables */
            napsolution (x, a, y, ay, ad, lo, hi, lambda, uf, kf, nuf, nkf,
                         loExists, hiExists, d_is_one) ;
        }
        Stat->nkf = nkf ;
    }

    /* wrapup and return solution */
    return (nap_wrapup (status, napdata, lambda)) ;
}

/* =========================================================================
 * ============================= napheap_setup =============================
 * ========================================================================= */
NAPdata * napheap_setup (void)
{
    int status ;
    NAPdata *Data ;
    status = NAPHEAP_STATUS_OK ;
    Data       = (NAPdata *) napheap_malloc (&status, 1, sizeof (NAPdata)) ;
    Data->Stat = (NAPstat *) napheap_malloc (&status, 1, sizeof (NAPstat)) ;
    Data->Parm = (NAPparm *) napheap_malloc (&status, 1, sizeof (NAPparm)) ;
    if ( status ) /* malloc failed */
    {
        Data = NULL ;
        return (Data) ;
    }
    napheap_default (Data->Parm) ;

    /* A starting guess for x is not required.  If the user does not
       provide a pointer to x, then NAPHEAP allocates memory for the
       solution x. napheap_terminate will free any allocated memory.
       The x_created pointer is only used internally by NAPHEAP when
       it creates memory for x */
    Data->x = NULL ;

    /* lambda, the multiplier associated with the constraint
       blo <= a'*x <= bhi, is initialized to +infinity. If the
       user gives a value to Data->lambda, then it is used as the
       starting guess for the multiplier. Otherwise, if Data->lambda
       remains +infinity, it is ignored. At the conclusion of the
       run, Data->lambda is the final multiplier associated with the
       linear constraint. */
    Data->lambda = NAPINF ;

    Data->n      = VSEMPTY ;/* problem dimension */
    Data->c      = NULL ;   /* linear term in the objective is -c'x */
    Data->d      = NULL ;   /* diagonal of objective Hessian */
    Data->a      = NULL ;   /* linear constraint vector */
    Data->blo    = -NAPINF; /* lower bound for a'x */
    Data->bhi    = NAPINF ; /* upper bound for a'x */
    Data->lo     = NULL ;   /* lower bound for x */
    Data->hi     = NULL ;   /* upper bound for x */
    Data->xWork  = NULL ;   /* NULL => let code allocate NAPFLOAT memory.
                                  Otherwise use xWork for NAPFLOAT memory.
                                  The amount of memory needed depends on the
                                  parameters in the Parm structure. See the
                                  memory allocation section of the code for the
                                  precise memory requirement (max size 5n). */
    Data->iWork  = NULL ;   /* NULL => let code allocate NAPINT memory.
                                  Otherwise use iWork for memory (4n+2) */

    /* the following are used internally when the code allocates memory */
    Data->x_created      = NULL ;
    Data->c_created      = NULL ;
    Data->xWork_created  = NULL ;
    Data->iWork_created  = NULL ;
    /* sum ak*ak/dk, k = 1:n, from a previous run */
    Data->akakdk         = NAPINF ;

    return (Data) ;
}

/* =========================================================================
 * ============================= napheap_terminate =========================
 * ========================================================================= */
void napheap_terminate
(
    NAPdata **DataHandle
)
{
    NAPdata *Data ;
    Data = *DataHandle ;
    if ( Data->x_created      != NULL ) napheap_free (Data->x_created) ;
    if ( Data->c_created      != NULL ) napheap_free (Data->c_created) ;
    if ( Data->xWork_created  != NULL ) napheap_free (Data->xWork_created) ;
    if ( Data->iWork_created  != NULL ) napheap_free (Data->iWork_created) ;
    napheap_free (Data->Stat) ;
    napheap_free (Data->Parm) ;
    napheap_free (Data) ;
}

/* ========================================================================== */
/* === non-user-callable functions ========================================== */
/* ========================================================================== */

/* The remainder of the functions below are not user-callable. */

/* ========================================================================== */
/* === nap_wrapup =========================================================== */
/* ========================================================================== */

/* finalize the solution, free workspace, and return to the user */

PRIVATE int nap_wrapup
(
    int         status,
    NAPdata   *napdata,
    NAPFLOAT    lambda
)
{
    NAPparm *Parm ;
    Parm = napdata->Parm ;
    napdata->lambda = lambda ;
    napdata->Stat->status = status ;
    if ( Parm->PrintStatus == TRUE )
    {
        napheap_print_status (napdata) ;
    }
    if ( Parm->PrintStat )
    {
        napheap_print_stat (napdata) ;
    }

    if ( !Parm->return_data )
    {
        if ( napdata->c_created != NULL )
        {
            if ( napdata->c_created == napdata->c )
            {
                napheap_free (napdata->c_created) ;
                napdata->c_created = napdata->c = NULL ;
            }
            else
            {
                napheap_free (napdata->c_created) ;
                napdata->c_created = NULL ;
            }
        }
        if ( napdata->xWork_created != NULL )
        {
            if ( napdata->xWork_created == napdata->xWork )
            {
                napheap_free (napdata->xWork_created) ;
                napdata->xWork_created = napdata->xWork = NULL ;
            }
            else
            {
                napheap_free (napdata->xWork_created) ;
                napdata->xWork_created = NULL ;
            }
        }
        if ( napdata->iWork_created != NULL )
        {
            if ( napdata->iWork_created == napdata->iWork )
            {
                napheap_free (napdata->iWork_created) ;
                napdata->iWork_created = napdata->iWork = NULL ;
            }
            else
            {
                napheap_free (napdata->iWork_created) ;
                napdata->iWork_created = NULL ;
            }
        }
    }
    return (status) ;
}

/*  ==========================================================================
    === napsearch ============================================================
    ==========================================================================
    If the slope of the dual function at lambda is positive, then search
    over break points to the right. If the slope is negative, then search
    over break points to the left.
    ========================================================================= */

PRIVATE int napsearch   /* return status of solution process */
(
    NAPFLOAT        *x, /* size n. solution (output) */
    NAPFLOAT   *lambda, /* optimal multiplier (input = guess, output = sol)*/
    NAPFLOAT   lambdaL, /* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT   lambdaR, /* optimal multiplier lies on [lambdaL, lambdaR] */
    NAPFLOAT         B, /* a'x = B */
    NAPFLOAT     slope, /* derivative of dual function */
    NAPFLOAT const *ay, /* size n. linear term in objective function yj/aj */
    NAPFLOAT const *ad, /* size n, diagonal of objective Hessian aj/dj */
    NAPFLOAT const  *a, /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPINT   const   n, /* problem dimension */
    NAPFLOAT *breakpts, /* size n, breakpts [j] is br_hi[j] or br_lo[j] */
    NAPFLOAT *breaknext,/* size n, TRUE if 2 break points in search dir. */
    NAPINT         nuf, /* number of unfixed variables */
    NAPINT *known_free, /* location of first known free variable in uf */
    NAPINT     *n_free, /* number of indices in free heap */
    NAPINT    *n_bound, /* number of indices in bound heap */
    NAPINT         *uf, /* size n, list of unfixed variables */
    NAPINT         *kf, /* size n, list of known free variables */
    NAPINT  *free_heap, /* size n+1, heap of free variables */
    NAPINT *bound_heap, /* size n+1, heap of bound variables */
    NAPstat      *Stat, /* solution statistics */
    NAPparm      *Parm, /* parameters */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists  /* TRUE when there are upper bound for x */
)
{
    NAPINT   jb, jfree, k, nb, nf, nbound, nfree, nkf, nk ;
    NAPFLOAT b, maxbound, maxfree, minbound,
             minfree, new_break, a2sum, a2sum_free, a2max ;

    /* initialization */
    b = B ;
    a2sum = NAPZERO ;      /* sum aj*aj/dj for all free variables */
    a2sum_free = NAPZERO ; /* sum aj*aj/dj for all known free variables */
    nbound = 0 ;           /* current number of bound variables */
    nfree = 0 ;            /* current number of free variables */
    nb = n ;               /* pointer to start of previously bound variables */
    nf = n ;               /* pointer to start of previously free variables */
    nkf = 0 ;              /* number of known free variables */
    const int d_is_one = Parm->d_is_one ; /* TRUE if d = 1 */
    const int a_is_one = (a == NULL) ? TRUE : FALSE ;
    if ( slope > NAPZERO ) /* move right */
    {
        /* d is positive and lambda is to the left of the true multiplier.
           Get ready to search over the break points to the right by
           constructing the heaps for currently bound and free variables,
           and by fixing variables whose break points are to the left.
           If the current iterate crossed the root, can exploit this in
           freeing variables. */
        LOOP_OVER_UNFIXED ;
        {
            const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
            if ( x [j] == VSZERO )
            {
                continue ;
            }
            else /* either 0 or 1 or 2 break points ahead */
            {
                const NAPFLOAT ayj = ay [j] ;
                const NAPFLOAT adj = aj ;
                const NAPFLOAT s = ayj - VSONE/adj ; /*lo break point*/
                if ( *lambda < s )            /* 2 break points ahead */
                {
                    nbound++ ;
                    bound_heap [nbound] = j ;
                    breakpts [j] = s ;
                    const NAPFLOAT r = ayj ; /*next brk */
                    if ( r < lambdaR )
                    {
                        breaknext [j] = r ; /* next break */
                    }
                    else
                    {
                        breaknext [j] = NAPINF ; /*no more brk pts*/
                    }
                }
                else /* 0 or 1 break points ahead */
                {
                    const NAPFLOAT r = ayj ; /* hi brk point*/
                    if ( r >= lambdaR )            /* j is known free */
                    {
                        a2sum_free += aj*adj ;
                        kf [nkf] = j ;
                        nkf++ ;
                        continue ;
                    }
                    else /* j free, exactly 1 break point ahead */
                    {
                        a2sum += aj*adj ;
                        nfree++ ;
                        free_heap [nfree] = j ;
                        breakpts [j] = r ;
                        breaknext [j] = NAPINF ; /* no more break pts */
                    }
                }
            }
        }
        END_UNFIXED_LOOP ;
        if ( (nuf == 0) && (nkf == 0) ) /* all variables bounds */
        {
            return (NAPHEAP_STATUS_INFEASIBLE) ;
        }

       /* recompute a2sum when it becomes much smaller than its max */
        a2sum += a2sum_free ;
        a2max = a2sum ;

        /* build heap for free variables */
        Stat->nfree = nfree ;
        napminheap_build (free_heap, breakpts, nfree) ;

        /* build heap for bound variables */
        Stat->nbound = nbound ;
        napminheap_build (bound_heap, breakpts, nbound) ;

        /* d is positive, lambda is to the left of the true multiplier, and
           the break points are stored in a heap.  Search over the break points
           to the right until finding 0 slope */
        if ( nfree > 0 )
        {
            jfree = free_heap [1] ;
            minfree = breakpts [jfree] ;
        }
        else minfree = NAPINF ;

        if ( nbound > 0 )
        {
            jb = bound_heap [1] ;
            minbound = breakpts [jb] ;
        }
        else minbound = NAPINF ;

        while ( 1 )
        {
            new_break = NAPMIN (minfree, minbound) ;
            const NAPFLOAT t = slope - a2sum*(new_break - *lambda) ;
            if ( (t <= NAPZERO) || (new_break >= NAPINF) ) /* done */
            {
                Stat->nbrks = (n-nf) + (n-nb) ;

                /* Final update when there are free variables */
                if ( (nfree > 0) || (nkf > 0) )
                {
                    *lambda += slope/a2sum ;
                    break ;
                }
                else if ( slope > NAPZERO ) /* Problem infeasible */
                {
                    return (NAPHEAP_STATUS_INFEASIBLE) ;
                }
            }
            *lambda = new_break ;
            slope = t ;

            /* update the heaps */
            while ( *lambda == minfree )
            {
                /* delete break point from heap */
                nfree = napminheap_delete (free_heap, breakpts, nfree) ;

                /* add break point to list of previously free variables */
                free_heap [nf] = jfree ;
                nf-- ;

                /* update a2sum */
                const NAPFLOAT aj = (a_is_one) ? VSONE : a [jfree] ;
                a2sum -= aj*aj ;

                /* store x [jfree], it must be bound */
                if ( aj > 0 ) x [jfree] = VSZERO ;
                else          x [jfree] = VSONE ;

                /* recompute a2sum if it has significant decay */
                if ( a2sum < Parm->decay*a2max )
                {
                    a2sum = a2sum_free ;
                    for (k = 1; k <= nfree; k++)
                    {
                        const NAPINT j = free_heap [k] ;
                        a2sum += (a_is_one) ? VSONE : a [j]*a [j] ;
                    }
                    a2max = a2sum ;
                }

                /* get the smallest entry in heap */
                if ( nfree == 0 )
                {
                    minfree = NAPINF ;
                    break ;
                }
                else
                {
                    jfree = free_heap [1] ;
                    minfree = breakpts [jfree] ;
                }
            }

            while ( *lambda == minbound )
            {
                /* delete break point from heap */
                nbound  = napminheap_delete(bound_heap, breakpts, nbound);

                /* add break point to list of previously bound  variables */
                bound_heap [nb] = jb ;
                nb-- ;

                /* compute new hi break point if necessary */
                const NAPFLOAT s = breaknext [jb] ;
                if ( s < NAPINF )
                {
                    breakpts [jb] = s ;

                    /* add new break point to free heap */
                    nfree = napminheap_add (free_heap, breakpts, jb, nfree);

                    /* get the smallest entry in free heap */
                    jfree = free_heap [1] ;
                    minfree = breakpts [jfree] ;
                    a2sum += (a_is_one) ? VSONE : a [jb]*a [jb] ;
                }
                else /* variable is known free */
                {
                    const NAPFLOAT r = (a_is_one) ? VSONE : a [jb]*a [jb] ;
                    a2sum_free += r ;
                    a2sum += r ;
                    kf [nkf] = jb ;
                    nkf++ ;
                }

                if ( a2sum > a2max )
                {
                    a2max = a2sum ;
                }

                /* get the smallest entry in bound heap */
                if ( nbound == 0 ) minbound = NAPINF ;
                else
                {
                    jb = bound_heap [1] ;
                    minbound = breakpts [jb] ;
                }
            }
        }
    }
    else if ( slope < NAPZERO ) /* move left */
    {
        /* d is positive and lambda is to the right of the true multiplier.
           Get ready to search over the break points to the left by
           constructing the heaps for currently bound and free variables,
           and by fixing variables whose break points are to the right.
           If the current iterate crossed the root, can exploit this to
           free variables. */
        LOOP_OVER_UNFIXED ;
        {
            const NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
            if ( x [j] == VSONE )
            {
                b -= aj ; /* fix j */
                continue ;
            }
            else /* either 0 or 1 or 2 break points ahead */
            {
                const NAPFLOAT ayj = ay [j] ;
                const NAPFLOAT adj = (a_is_one) ? VSONE : aj ;
                const NAPFLOAT s = ayj ; /* hi brk point */
                if ( *lambda > s )        /* 2 break points ahead */
                {
                    nbound++ ;
                    bound_heap [nbound] = j ;
                    breakpts [j] = s ;
                    const NAPFLOAT r = ayj - VSONE/adj ;/*next brk */
                    if ( r > lambdaL )
                    {
                            breaknext [j] = r ; /* next break */
                    }
                    else
                    {
                        breaknext [j] = NAPINF ; /*no more brk pts*/
                    }
                }
                else /* 0 or 1 break points ahead */
                {
                    const NAPFLOAT r = ayj - VSONE/adj ; /* lo brk point*/
                    if ( r <= lambdaL ) /* j is known free */
                    {
                        a2sum_free += aj*adj ;
                        kf [nkf] = j ;
                        nkf++ ;
                        continue ;
                    }
                    else /* j free, exactly 1 break point ahead */
                    {
                        a2sum += aj*adj ;
                        nfree++ ;
                        free_heap [nfree] = j ;
                        breakpts [j] = r ;
                        breaknext [j] = NAPINF ; /* no more break pts */
                    }
                }
            }
        }
        END_UNFIXED_LOOP ;
        if ( (nuf == 0) && (nkf == 0) ) /* all variables at bounds */
        {
            return (NAPHEAP_STATUS_INFEASIBLE) ;
        }

       /* recompute a2sum when it becomes much smaller than its max */
        a2sum += a2sum_free ;
        a2max = a2sum ;

        /* build the free heap */
        Stat->nfree = nfree ;
        napmaxheap_build (free_heap, breakpts, nfree) ;

        /* build the bound heap */
        Stat->nbound = nbound ;
        napmaxheap_build (bound_heap, breakpts, nbound) ;

        /* d is positive, lambda is to the right of the true multiplier, and
           the break points are stored in a heap.  Search over the break points
           to the left until finding 0 slope */
        if ( nfree > 0 )
        {
            jfree = free_heap [1] ;
            maxfree = breakpts [jfree] ;
        }
        else maxfree = -NAPINF ;

        if ( nbound > 0 )
        {
            jb = bound_heap [1] ;
            maxbound = breakpts [jb] ;
        }
        else maxbound = -NAPINF ;

        while ( 1 )
        {
            new_break = NAPMAX (maxfree, maxbound) ;
            const NAPFLOAT t = slope - a2sum*(new_break - *lambda) ;
            if ( (t >= NAPZERO) || (new_break <= -NAPINF) ) /* done */
            {
                Stat->nbrks = (n-nb) + (n-nf) ;

                /* Final update when there are free variables */
                if ( ((nfree > 0) || (nkf > 0)) && (a2sum > NAPZERO) )
                {
                    *lambda += slope/a2sum ;
                    break ;
                }
                else if ( slope < NAPZERO ) /* problem is infeasible */
                {
                    return (NAPHEAP_STATUS_INFEASIBLE) ;
                }
            }
            *lambda = new_break ;
            slope = t ;

            /* update the heaps */
            while ( *lambda == maxfree )
            {
                /* delete break point from heap */
                nfree = napmaxheap_delete (free_heap, breakpts, nfree) ;

                /* add break point to list of previously free variables */
                free_heap [nf] = jfree ;
                nf-- ;

                /* update a2sum */
                const NAPFLOAT aj = (a_is_one) ? VSONE : a [jfree] ;
                a2sum -= (a_is_one) ? VSONE : aj*aj ;

                /* recompute a2sum if it has significant decay */
                if ( a2sum < Parm->decay*a2max )
                {
                    a2sum = a2sum_free ;
                    for (k = 1; k <= nfree; k++)
                    {
                        const NAPINT j = free_heap [k] ;
                        a2sum += (a_is_one) ? VSONE : a [j]*a [j] ;
                    }
                    a2max = a2sum ;
                }

                /* store x [jfree] */
                x [jfree] = VSONE ;

                /* get the smallest entry in heap */
                if ( nfree == 0 )
                {
                    maxfree = -NAPINF ;
                    break ;
                }
                else
                {
                    jfree = free_heap [1] ;
                    maxfree = breakpts [jfree] ;
                }
            }
            while ( *lambda == maxbound )
            {
                /* delete break point from heap */
                nbound  = napmaxheap_delete(bound_heap, breakpts, nbound);

                /* add break point to list of previously bound variables */
                bound_heap [nb] = jb ;
                nb-- ;

                /* compute new lo break point if necessary */
                const NAPFLOAT s = breaknext [jb] ;
                if ( s < NAPINF )
                {
                    breakpts [jb] = s ;

                    /* add new break point to free heap */
                    nfree = napmaxheap_add (free_heap, breakpts, jb, nfree);

                    /* get the smallest entry in free heap */
                    jfree = free_heap [1] ;
                    maxfree = breakpts [jfree] ;
                    a2sum += (a_is_one) ? VSONE : a [jb]*a [jb] ;
                }
                else /* variable is known free */
                {
                    const NAPFLOAT r = (a_is_one) ? VSONE : a [jb]*a [jb] ;
                    a2sum_free += r ;
                    a2sum += r ;
                    kf [nkf] = jb ;
                    nkf++ ;
                }

                if ( a2sum > a2max )
                {
                    a2max = a2sum ;
                }

                /* get the largest entry in bound heap */
                if ( nbound == 0 ) maxbound = -NAPINF ;
                else
                {
                    jb = bound_heap [1] ;
                    maxbound = breakpts [jb] ;
                }
            }
        }
    }

    *n_free  = nfree;
    *n_bound = nbound;
    *known_free = nkf ;
    return (NAPHEAP_STATUS_OK) ;
}

/*  =========================================================================
    === napsolution =========================================================
    =========================================================================
    d is positive and the optimal dual solution lambda has been determined.
    Evaluate the optimal primal solution.
    =========================================================================*/

PRIVATE void napsolution
(
    NAPFLOAT        *x, /* size n, solution (output) */
    NAPFLOAT const  *a, /* size n, coefficient of linear constraint */
    NAPFLOAT const  *y, /* size n, linear term in objective */
    NAPFLOAT const *ay, /* size n, linear term y in objective over a */
    NAPFLOAT const *ad, /* size n, coefficient of linear constraint over d */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPINT   const *uf, /* size n, unfixed variables */
    NAPINT   const *kf, /* size n, known free variables */
    NAPINT         nuf, /* number of unfixed variables */
    NAPINT         nkf, /* number of known free variables */
    int const loExists, /* TRUE when there are lower bound for x */
    int const hiExists, /* TRUE when there are upper bound for x */
    int const d_is_one  /* TRUE when diagonal is identically one */
)
{
    NAPFLOAT t ;
    NAPINT k ;

    int const a_is_one = (a == NULL) ? TRUE : FALSE ;
    /* set unbound variables */
    LOOP_OVER_UNFIXED_NO_UF_UPDATE ;
        t = (a_is_one) ? y [j] - lambda : y [j] - lambda*a [j] ;
        if      ( t < VSZERO ) t = VSZERO ;
        else if ( t > VSONE )  t = VSONE ;
        x [j] = t ;
    END_LOOP
    
    /* set free variables */
    for (k = 0; k < nkf; k++)
    {
        NAPINT j = kf [k] ;
        x [j] = (a_is_one) ? y [j] - lambda : y [j] - lambda*a [j];
    }
}

/*  =========================================================================
    === napsolution1 ========================================================
    =========================================================================
    d is positive and the optimal dual solution lambda has been determined.
    Evaluate the optimal primal solution. Sets variables in free heap and
    bound heap to optimal values. Assumes napsearch has been called.
    =========================================================================*/

PRIVATE void napsolution1
(
    NAPFLOAT        *x,  /* size n. solution (output) */
    NAPFLOAT const  *y,  /* size n. linear term in objective function */
    NAPFLOAT const  *d,  /* size n, diagonal of objective Hessian */
    NAPFLOAT const  *a,  /* size n, linear constraint vector */
    NAPFLOAT const *lo, /* size n, lower bounds on x */
    NAPFLOAT const *hi, /* size n, upper bounds on x */
    NAPFLOAT    lambda, /* optimal multiplier */
    NAPFLOAT     slope, /* slope of dual function determines search direction */
    NAPINT   const *kf,/* size n, known free variables */
    NAPINT         nkf, /* number of known free variables */
    NAPINT  *free_heap,
    NAPINT *bound_heap,
    NAPINT      n_free,
    NAPINT     n_bound,
    int const d_is_one
)
{
    NAPINT i ;
    NAPFLOAT const Lambda = lambda ;

    int const a_is_one = (a == NULL) ? TRUE : FALSE ;
    /* set variables in bound heap */
    for (i = 1; i <= n_bound; i++)
    {
        NAPINT j = bound_heap [i] ;
        if ( slope < NAPZERO ) x [j] = VSZERO ; /* searching left */
        else                   x [j] = VSONE ;
    }

    /* set variables in free heap */
    for (i = 1; i <= n_free; i++)
    {
        NAPINT j = free_heap [i] ;
        x [j] = (a_is_one) ? y [j] - Lambda : y [j] - Lambda*a [j] ;
    }
    /* set additional free variables */
    for (i = 0; i < nkf; i++)
    {
        NAPINT j = kf [i] ;
        x [j] = (a_is_one) ? y [j] - Lambda : y [j] - Lambda*a [j] ;
    }
}

/* ========================================================================== */
/* === nap_work ============================================================= */
/* ========================================================================== */

/* allocate NAPFLOAT workspace */

PRIVATE int nap_work
(
    NAPINT  const      n,
    NAPFLOAT      *xWork,     /* user-provide workspace */
    NAPFLOAT       **xw2,     /* return malloc'ed space */

    /* the 6 work arrays, each of size n */
    NAPFLOAT        **ad,
    NAPFLOAT        **ay,
    NAPFLOAT     **br_hi,
    NAPFLOAT     **br_lo,
    NAPFLOAT  **breakpts,
    NAPFLOAT **breaknext
)
{
    NAPFLOAT *work ;

    /* determine how many workspaces to allocate */
    int m, status ;

    status = NAPHEAP_STATUS_OK ;
    if ( xWork == NULL )
    {
        m = (ad != NULL) + (ay != NULL) + (br_hi != NULL) +
            (br_hi != NULL) + (breakpts != NULL) + (breaknext != NULL) ;
        /* n will not be zero at this point */
        (*xw2) = (NAPFLOAT *) napheap_malloc (&status, NAPMAX(m*n,1),
                                              sizeof (NAPFLOAT)) ;
        if ( status )
        {
            return (NAPHEAP_STATUS_OUT_OF_MEMORY) ;
        }
        work = (*xw2) ;
    }
    else
    {
        /* use the user-provided workspace */
        work = xWork ;
        (*xw2) = NULL ;
    }

    /* allocate each array.  Each array is size n */
    if (breaknext != NULL)
    {
        (*breaknext) = work ;   work += n ;
    }
    if (breakpts != NULL)
    {
        (*breakpts) = work ;    work += n ;
    }
    if (ad != NULL)
    {
        (*ad) = work ;          work += n ;
    }
    if (ay != NULL)
    {
        (*ay) = work ;          work += n ;
    }
    if (br_hi != NULL)
    {
        (*br_hi) = work ;       work += n ;
    }
    if (br_lo != NULL)
    {
        (*br_lo) = work ;       work += n ;
    }

    return (status) ;
}

/* ========================================================================== */
/* === napheap_default ====================================================== */
/* ========================================================================== */

/*  Set default parameter values for the napheap code. If the parameter
    argument of napheap is NULL, this routine is called automatically to
    set parameters to their default value.  If the user wishes to set
    parameter values, then the NAPparm structure should be declared
    in the main program. The user could call napheap_default to initialize
    the structure, and then individual elements in the structure
    could be changed, before passing the structure to napheap.
*/

void napheap_default
(
    NAPparm *Parm
)
{
    /* ---------------------------------------------------------------------- */
    /* Parm->PrintStatus */
    /* T => Print the status of the run */
    Parm->PrintStatus = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->PrintStat */
    /* T => Print the statistics for the run */
    Parm->PrintStat = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->PrintParm */
    /* T => Print the parameter values used for the run */
    Parm->PrintParm = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->use_prior_data */
    /* T => Use malloc'd arrays and data from a previous run.
            These arrays were saved during prior run by setting
            parameter return_data = TRUE */
    Parm->use_prior_data = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->return_data */
    /* T => Return malloc'd arrays for the current run in the
            input data structure napdata. This is useful when solving
            another problem where only the linear input data y changes. */
    Parm->return_data = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->loExists */
    /* T => the user provides lower bounds for the variables
       F => treat the lower bounds as 0
            Note: The lo argument of napheap is ignored, it could be NULL.
                  If the lo argument is NULL, then loExists is treated as F */
    Parm->loExists = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->hiExists */
    /* T => the user provides upper bounds for the variables
       F => treat the upper bounds as 1
            Note: The hi argument of napheap is ignored, it could be NULL.
                  If the hi argument is NULL, then hiExists is treated as F */
    Parm->hiExists = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->Aexists */
    /* T => there is a linear constraint bl <= a'x <= bu
       F => treat the upper bounds as +infinity */
    Parm->Aexists = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->d_is_pos: */
    /* T => all entries of diagonal of objective Hessian are > 0
       F => there could be zeros on the diagonal of the Hessian. */
    Parm->d_is_pos = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->d_is_zero: */
    /* T => all entries of diagonal of objective Hessian are zero
            Note: The d argument of napheap is ignored. It could be NULL
       F => there could be nonzeros on the diagonal of the Hessian. */
    Parm->d_is_zero = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->d_is_one: */
    /* T => diagonal of the objective Hessian is identically 1
            Note: The d argument of napheap is ignored. It could be NULL */
    Parm->d_is_one = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->K: */
    /* Upper bound on the number of Newton or variable fixing iterations.
       After performing these K iterations, the code switches to the break
       point search. The Newton algorithm also switches to the breakpoint
       search if two iterates are generated on opposite sides of the root. */
    Parm->K = 20 ;

    /* ---------------------------------------------------------------------- */
    /* Parm->newton: */
    /* T => use the Newton to update lambda, switch to the break point
            search either at iteration K or at the first time that two
            Newton iterates lie on opposite side of the root (whichever
            happens first)
       F => use K iterations of the variable fixing algorithm before
            switching to a breakpoint search */
    Parm->newton = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->newton_scale: */
    /* An ordinary Newton iteration often converges monotonically to the
       root. Hence, the Newton iterates are scaled by a factor greater than 1
       in an effort to produce an iterate on the opposite side of the root.
       When such an iterate is generated, the algorithm switches to a
       breakpoint search. The scaled iteration is given by

       lambda_new = lambda_old - newton_scale*L'(lambda_old)/L''(lambda_old) */

    Parm->newton_scale = 1.10 ;

    /* ---------------------------------------------------------------------- */
    /* Parm->decay: */
    /* To combat the effect of rounding errors, a2sum is recomputed whenever
       a2sum <= decay*a2max, the maximum a2sum encountered */
    Parm->decay = 0.01 ;

    /* ---------------------------------------------------------------------- */
    /* Parm->refine: */
    /* T => perform an iterative refinement step to improve solution accuracy
       F => no refinement
       Since the accuracy has been improved by monitoring the size of the
       a2sum variable, this parameter is ignored when d is positive.  In the
       next version of the code, this parameter will be removed. */
    Parm->refine = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* Parm->err: */
    /* The solution is refined if |a'x - b| > err. This parameter will also
       be removed in the next version of the code. */
    Parm->err = NAPZERO ;

    /* ---------------------------------------------------------------------- */
    /* Parm->check: */
    /* Check the input parameters for dj < 0 or hij < loj or Bhi < Blo */
    Parm->check = FALSE ;
}

/* ========================================================================== */
/* === Source/napheap_print.c =============================================== */
/* ========================================================================== */

/* user-callable functions for printing napheap parameters and statistics. */

#ifdef MATLAB_MEX_FILE
#define LPAREN "("
#define RPAREN ")"
#define OFFSET 1
#else
#define LPAREN "["
#define RPAREN "]"
#define OFFSET 0
#endif

/* ========================================================================== */
/* === napheap_print_status ================================================= */
/* ========================================================================== */

/* print the status at termination of the napheap run */

void napheap_print_status
(
    NAPdata *Data
)
{
    NAPstat *Stat = Data->Stat ;
    int status = Stat->status ;
    printf ("\nNAPHEAP run status (Version %d.%d.%d, %s):\n\n",
        NAPHEAP_MAIN_VERSION, NAPHEAP_SUB_VERSION, NAPHEAP_SUBSUB_VERSION,
        NAPHEAP_DATE) ;

    if ( status == NAPHEAP_STATUS_OK )
    {
        printf ("No errors detected by NAPHEAP when solving "
                "the knapsack problem.\n") ;
    }
    else if ( status == NAPHEAP_STATUS_INFEASIBLE )
    {
        printf ("Knapsack problem in NAPHEAP is infeasible.\n") ;
    }
    else if ( status == NAPHEAP_STATUS_UNBOUNDED )
    {
        printf ("Optimal value for the knapsack problem in NAPHEAP "
                "is -infinity.\n") ;
    }
    else if ( status == NAPHEAP_STATUS_OUT_OF_MEMORY )
    {
        printf ("NAPHEAP ran out of memory.\n") ;
    }
    else if ( status == NAPHEAP_STATUS_INVALID_N )
    {
        printf ("In NAPHEAP, the specified problem dimension "
                "is either missing or nonpositive (%ld).\n",
                (int64_t) Stat->kerror) ;
    }
    else if ( status == NAPHEAP_STATUS_INVALID_BOUNDS )
    {
        printf ("In NAPHEAP, the problem bounds are invalid since "
                "lo %s%ld%s = %e > hi %s%ld%s = %e.\n",
            LPAREN, (int64_t) Stat->kerror + OFFSET, RPAREN, Stat->lobad,
            LPAREN, (int64_t) Stat->kerror + OFFSET, RPAREN, Stat->hibad) ;
    }
    else if ( status == NAPHEAP_STATUS_INVALID_D )
    {
        if ( Stat->dbad < 0 )
        {
            printf ("In NAPHEAP, the Hessian diagonal d is invalid "
                    "(d %s%ld%s = %e < 0).\n", LPAREN,
                     (int64_t) Stat->kerror + OFFSET, RPAREN, Stat->dbad) ;
        }
        else
        {
            printf ("In NAPHEAP, the parameter Parm.d_is_pos is TRUE "
                    "but d %s%ld%s = 0\n", LPAREN,
                    (int64_t) Stat->kerror + OFFSET, RPAREN) ;
        }
    }
    else if ( status ==  NAPHEAP_STATUS_CONFLICTING_PARM_D )
    {
        printf ("The parameters d_is_one, d_is_pos, and d_is_zero are "
                "assigned conflicting values.\n") ;
        printf ("For example, both d_is_zero and d_is_one are TRUE, which is "
                "impossible.\n") ;
        printf ("If d_is_one and d_is_zero are both FALSE, then the d argument "
                "of napheap cannot be NULL.\n") ;
    }
    else if ( status == NAPHEAP_STATUS_MISSING_PARM )
    {
        printf ("In NAPHEAP, the Parm structure of the input napdata\n"
                "is missing. This should be created by napheap_setup.") ;
    }
    else if ( status == NAPHEAP_STATUS_MISSING_STAT )
    {
        printf ("In NAPHEAP, the Stat structure of the input napdata\n"
                "is missing. This should be created by napheap_setup.") ;
    }
}

/* ========================================================================== */
/* === napheap_print_stat =================================================== */
/* ========================================================================== */

/* print the napheap statistics */
void napheap_print_stat
(
    NAPdata *Data
)
{
    NAPstat *Stat = Data->Stat ;
    if (Stat == NULL)
    {
        return ;
    }

    printf ("\nNAPHEAP statistics (Version %d.%d.%d, %s):\n\n",
        NAPHEAP_MAIN_VERSION, NAPHEAP_SUB_VERSION, NAPHEAP_SUBSUB_VERSION,
        NAPHEAP_DATE) ;

    printf ("Stat.nkf (number of known free variables):           %ld\n",
        (int64_t) Stat->nkf) ;
    printf ("Stat.nfree (free variables in initial heap):         %ld\n",
        (int64_t) Stat->nfree) ;
    printf ("Stat.nbound (bound variables in initial heap):       %ld\n",
        (int64_t) Stat->nbound) ;
    printf ("Stat.nbrks (break points to reach initial solution): %ld\n",
        (int64_t) Stat->nbrks) ;
    printf ("Stat.nrefine (break points during refinement):       %ld\n",
        (int64_t) Stat->nrefine) ;
    printf ("Stat.nvarfix (variable fixing iterations):           %ld\n",
        (int64_t) Stat->nvarfix) ;
    printf ("Stat.nnewton (number of Newton iterations):          %ld\n",
        (int64_t) Stat->nnewton) ;
    printf ("Stat.nsecant (number of secant iterations):          %ld\n",
        (int64_t) Stat->nsecant) ;
    printf ("\n") ;
}

/* ========================================================================== */
/* === napheap_print_parm ================================================== */
/* ========================================================================== */

void napheap_print_parm
(
    NAPdata  *Data        /* parameter structure */
)
{
    NAPparm *Parm, Default_Parm ;

    if ( Data == NULL ) Parm = NULL ;
    else                Parm = Data->Parm ;
    printf ("\nNAPHEAP parameter settings (Version %d.%d.%d, %s):\n",
        NAPHEAP_MAIN_VERSION, NAPHEAP_SUB_VERSION, NAPHEAP_SUBSUB_VERSION,
        NAPHEAP_DATE) ;

    if (Parm == NULL)
    {
        /* initialize the default parameters */
        Parm = &Default_Parm ;
        napheap_default (Parm) ;
    }

    printf ("Parm.PrintStatus    ") ;
    if ( Parm->PrintStatus == TRUE )
    {
        printf ("TRUE (print the status of the run)\n") ;
    }
    else
    {
        printf ("FALSE (do not print the status of the run)\n") ;
    }

    printf ("Parm.PrintStat      ") ;
    if ( Parm->PrintStat == TRUE )
    {
        printf ("TRUE (print the statistics for the run)\n") ;
    }
    else
    {
        printf ("FALSE (do not print the statistics for the run)\n") ;
    }

    printf ("Parm.PrintParm      ") ;
    if ( Parm->PrintParm == TRUE )
    {
        printf ("TRUE (print the parameters used for the run)\n") ;
    }
    else
    {
        printf ("FALSE (do not print the parameters used for the run)\n") ;
    }

    printf ("Parm.use_prior_data ") ;
    if ( Parm->use_prior_data == TRUE )
    {
        printf ("TRUE (use internal arrays from previous run)\n") ;
    }
    else
    {
        printf ("FALSE (do not use arrays from previous run)\n") ;
    }

    printf ("Parm.return_data    ") ;
    if ( Parm->return_data == TRUE )
    {
        printf ("TRUE (return arrays from current run in napdata input\n"
                "                           argument of napheap)\n") ;
    }
    else
    {
        printf ("FALSE (do not return data arrays from current run)\n") ;
    }

    printf ("Parm.loExists       ") ;
    if ( Parm->loExists == TRUE )
    {
        printf ("TRUE (lower bounds on the variables are present)\n") ;
    }
    else
    {
        printf ("FALSE (lower bounds on the variables are treated as "
                "-infinity)\n") ;
    }

    printf ("Parm.hiExists       ") ;
    if ( Parm->hiExists == TRUE )
    {
        printf ("TRUE (upper bounds on the variables are present)\n") ;
    }
    else
    {
        printf ("FALSE (upper bounds on the variables are treated as "
                "+infinity)\n") ;
    }

    printf ("Parm.Aexists        ") ;
    if ( Parm->Aexists == TRUE )
    {
        printf ("TRUE (the linear constraint bl <= a'x <= bu is present)\n") ;
    }
    else
    {
        printf ("FALSE (the only constraints are bound constraints "
                "lo <= x <= hi)\n") ;
    }

    printf ("Parm.d_is_zero:     ") ;
    if ( Parm->d_is_zero )
    {
        printf ("TRUE (diagonal of objective Hessian is zero)\n") ;
        printf ("                          d argument of napheap and "
                                    "Parm.d_is_pos are ignored\n") ;
    }
    else
    {
        printf ("FALSE (d, diagonal of objective Hessian is nonzero)\n");
    }

    printf ("Parm.d_is_one:      ") ;
    if ( Parm->d_is_one )
    {
        printf ("TRUE (diagonal of objective Hessian is identically one)\n") ;
        printf ("                           d argument of napheap and "
                                    "Parm.d_is_pos are ignored\n") ;
    }
    else
    {
        printf ("FALSE (d, diagonal of objective Hessian is not all one)\n") ;
    }

    printf ("Parm.d_is_pos:      ") ;
    if ( Parm->d_is_pos )
    {
        printf ("TRUE (diagonal of objective Hessian is strictly positive)\n") ;
        printf ("                          Parm.d_is_pos ignored if either "
                                    "Parm.d_is_zero or\n") ;
        printf ("                          Parm.d_is_one are TRUE\n") ;
    }
    else
    {
        printf ("FALSE (d, diagonal of objective Hessian can have either\n") ;
        printf ("                                   "
                "positive or zero entries)\n") ;
    }

    printf ("Parm.refine:        ") ;
    if ( Parm->refine )
    {
        printf ("TRUE (perform refinement step)\n") ;
        printf ("Parm.err:           %e (solution refined if constraint"
            " violation exceeds this value)\n", Parm->err) ;
    }
    else
    {
        printf ("FALSE (do not perform refinement step)\n") ;
        printf ("Parm.err:           ignored (only used for refinement)\n") ;
    }

    printf ("Parm.K:             ") ;
    if ( Parm->K <= 0 )
    {
        printf ("0 (only use break point searching algorithm)\n") ;
        printf ("Parm.newton:        ignored (only used if K > 1)\n") ;
    }
    else
    {
        printf ("%d (max number of Newton or variable fixing iterations)\n",
                Parm->K) ;
        printf ("Parm.newton:        ") ;
        if ( Parm->newton )
        {
            printf ("TRUE (use Newton's method)\n") ;
            printf ("Parm.newton_scale:  ") ;
            printf ("%6.4f (multiply Newton step by this factor)\n",
                    Parm->newton_scale) ;
        }
        else
        {
            printf ("FALSE (use variable fixing algorithm)\n") ;
        }
    }

    printf ("Parm.decay:         ") ;
    printf ("%6.4f (recompute L'' when it decreases by decay factor)\n",
             Parm->decay) ;

    printf ("Parm.check:         ") ;
    if ( Parm->check )
    {
        printf ("TRUE (check the problem data for consistency)\n") ;
    }
    else
    {
        printf ("FALSE (do not check the problem data)\n") ;
    }
}

/* construct the minheap functions */
#define CONSTRUCT_MINHEAP_FUNCTIONS
/* ========================================================================== */
/* === napminmaxheap.c ====================================================== */
/* ========================================================================== */

/* This module creates two kinds of heap: a min heap and a maxheap */

#ifdef CONSTRUCT_MINHEAP_FUNCTIONS
/* min heap functions */
#define NAPHEAP_BUILD  napminheap_build
#define NAPHEAP_DELETE napminheap_delete
#define NAPHEAP_ADD    napminheap_add
#define NAPHEAPIFY     napminheapify
#define IN_ORDER(a,b)     (a < b)
#define OUT_OF_ORDER(a,b) (a > b)
#endif

#ifdef CONSTRUCT_MAXHEAP_FUNCTIONS
/* max heap functions */
#define NAPHEAP_BUILD  napmaxheap_build
#define NAPHEAP_DELETE napmaxheap_delete
#define NAPHEAP_ADD    napmaxheap_add
#define NAPHEAPIFY     napmaxheapify
#define IN_ORDER(a,b)     (a > b)
#define OUT_OF_ORDER(a,b) (a < b)
#endif

/* For the macros above:
    in_order(a,b) is true if [a,b] is in the proper order in the heap.
    out_of_order(a,b) is true if [a,b] is not in proper order in the heap */

/* ==========================================================================
   === NAPHEAP_BUILD ========================================================
   ==========================================================================
   build a min/max heap in heap [1..nheap]
   heap [i] is an index of an element of array X
   ========================================================================== */

PRIVATE void NAPHEAP_BUILD
(
    NAPINT      *heap, /* on input, unsorted set of indices.  size nheap+1 */
    NAPFLOAT const *X, /* numbers to sort (not modified) */
    NAPINT      nheap  /* number of elements to build into the heap */
)
{
    NAPINT p ;

    for (p = nheap/2 ; p >= 1 ; p--)
    {
        NAPHEAPIFY (heap, X, p, nheap) ;
    }
}

/* ==========================================================================
   === NAPHEAP_DELETE =======================================================
   ==========================================================================
   delete the top element in a min/max heap
   ========================================================================== */

PRIVATE NAPINT NAPHEAP_DELETE  /* return new size of heap */
(
    NAPINT      *heap, /* indices into X, 1..n on input.  size nheap+1 */
    NAPFLOAT const *X, /* not modified */
    NAPINT      nheap  /* number of items in heap */
)
{
    if (nheap <= 1)
    {
        return (0) ;
    }

    /* move element from the end of the heap to the top */
    heap [1] = heap [nheap] ;
    nheap-- ;
    NAPHEAPIFY (heap, X, 1, nheap) ;
    return (nheap) ;
}

/* ========================================================================== */
/* === NAPHEAP_ADD ========================================================== */
/* ========================================================================== */

/* add a new leaf to a min/max heap */

PRIVATE NAPINT NAPHEAP_ADD      /* returns the new size of the heap */
(
    NAPINT      *heap, /* size n, containing indices into X.  size nheap+2 */
    NAPFLOAT const *X, /* not modified */
    NAPINT       leaf, /* the new leaf */
    NAPINT      nheap  /* number of elements in heap not counting new one */
)
{
    NAPFLOAT Xold, Xnew ;
    NAPINT i, newindex, old ;

    nheap++ ;
    old = nheap ;
    heap [old] = leaf ;
    Xold = X [leaf] ;
    while ( old > 1 )
    {
        newindex = old/2 ;
        i = heap [newindex] ;
        Xnew = X [i] ;
        if (OUT_OF_ORDER (Xnew, Xold)) /* for min heap: ( Xnew > Xold ) */
        {
            /* Xnew and Xold are out of order, so swap them */
            heap [newindex] = leaf ;
            heap [old] = i ;
        }
        else return (nheap) ;
        old = newindex ;
    }
    return (nheap) ;
}

/* ========================================================================== */
/* === napminheapify ======================================================== */
/* ========================================================================== */

/* heapify starting at node p.  On input, the heap at node p satisfies the */
/* heap property, except for heap [p] itself.  On output, the whole heap */
/* satisfies the heap property. */

PRIVATE void NAPHEAPIFY
(
    NAPINT      *heap, /* size nheap+1, containing indices into X */
    NAPFLOAT const *X, /* not modified */
    NAPINT          p, /* start at node p in the heap */
    NAPINT      nheap  /* heap [1 ... n] is in use */
)
{
    NAPFLOAT Xe, Xleft, Xright ;
    NAPINT left, right, e, hleft, hright ;

    e = heap [p] ;
    Xe = X [e] ;

    while ( TRUE )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            /* both left and right children of p are in the heap */
            hleft  = heap [left] ;
            hright = heap [right] ;
            Xleft  = X [hleft] ;
            Xright = X [hright] ;
            if (IN_ORDER (Xleft, Xright)) /* for min heap: (Xleft < Xright) */
            {
                /* traverse down the left heap */
                if (OUT_OF_ORDER (Xe, Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
            else
            {
                /* traverse down the right heap */
                if (OUT_OF_ORDER (Xe, Xright)) /* for min heap: (Xe > Xright) */
                {
                    /* swap e with the right child */
                    heap [p] = hright ;
                    p = right ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
        }
        else
        {
            /* only the left child of p is in the heap */
            if ( left <= nheap )
            {
                /* traverse down the left heap */
                hleft = heap [left] ;
                Xleft = X [hleft] ;
                if (OUT_OF_ORDER (Xe, Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
            }
            /* element e, with value Xe, has found its final home */
            heap [p] = e ;
            return ;
        }
    }
}

#undef NAPHEAP_BUILD
#undef NAPHEAP_DELETE
#undef NAPHEAP_ADD
#undef NAPHEAPIFY
#undef IN_ORDER
#undef OUT_OF_ORDER
#undef CONSTRUCT_MINHEAP_FUNCTIONS
#undef CONSTRUCT_MAXHEAP_FUNCTIONS

/* construct the maxheap functions */
#define CONSTRUCT_MAXHEAP_FUNCTIONS

/* ========================================================================== */
/* === napminmaxheap.c ====================================================== */
/* ========================================================================== */

/* This module creates two kinds of heap: a min heap and a maxheap */

#ifdef CONSTRUCT_MINHEAP_FUNCTIONS
/* min heap functions */
#define NAPHEAP_BUILD  napminheap_build
#define NAPHEAP_DELETE napminheap_delete
#define NAPHEAP_ADD    napminheap_add
#define NAPHEAPIFY     napminheapify
#define IN_ORDER(a,b)     (a < b)
#define OUT_OF_ORDER(a,b) (a > b)
#endif

#ifdef CONSTRUCT_MAXHEAP_FUNCTIONS
/* max heap functions */
#define NAPHEAP_BUILD  napmaxheap_build
#define NAPHEAP_DELETE napmaxheap_delete
#define NAPHEAP_ADD    napmaxheap_add
#define NAPHEAPIFY     napmaxheapify
#define IN_ORDER(a,b)     (a > b)
#define OUT_OF_ORDER(a,b) (a < b)
#endif

/* For the macros above:
    in_order(a,b) is true if [a,b] is in the proper order in the heap.
    out_of_order(a,b) is true if [a,b] is not in proper order in the heap */

/* ==========================================================================
   === NAPHEAP_BUILD ========================================================
   ==========================================================================
   build a min/max heap in heap [1..nheap]
   heap [i] is an index of an element of array X
   ========================================================================== */

PRIVATE void NAPHEAP_BUILD
(
    NAPINT      *heap, /* on input, unsorted set of indices.  size nheap+1 */
    NAPFLOAT const *X, /* numbers to sort (not modified) */
    NAPINT      nheap  /* number of elements to build into the heap */
)
{
    NAPINT p ;

    for (p = nheap/2 ; p >= 1 ; p--)
    {
        NAPHEAPIFY (heap, X, p, nheap) ;
    }
}

/* ==========================================================================
   === NAPHEAP_DELETE =======================================================
   ==========================================================================
   delete the top element in a min/max heap
   ========================================================================== */

PRIVATE NAPINT NAPHEAP_DELETE  /* return new size of heap */
(
    NAPINT      *heap, /* indices into X, 1..n on input.  size nheap+1 */
    NAPFLOAT const *X, /* not modified */
    NAPINT      nheap  /* number of items in heap */
)
{
    if (nheap <= 1)
    {
        return (0) ;
    }

    /* move element from the end of the heap to the top */
    heap [1] = heap [nheap] ;
    nheap-- ;
    NAPHEAPIFY (heap, X, 1, nheap) ;
    return (nheap) ;
}

/* ========================================================================== */
/* === NAPHEAP_ADD ========================================================== */
/* ========================================================================== */

/* add a new leaf to a min/max heap */

PRIVATE NAPINT NAPHEAP_ADD      /* returns the new size of the heap */
(
    NAPINT      *heap, /* size n, containing indices into X.  size nheap+2 */
    NAPFLOAT const *X, /* not modified */
    NAPINT       leaf, /* the new leaf */
    NAPINT      nheap  /* number of elements in heap not counting new one */
)
{
    NAPFLOAT Xold, Xnew ;
    NAPINT i, newindex, old ;

    nheap++ ;
    old = nheap ;
    heap [old] = leaf ;
    Xold = X [leaf] ;
    while ( old > 1 )
    {
        newindex = old/2 ;
        i = heap [newindex] ;
        Xnew = X [i] ;
        if (OUT_OF_ORDER (Xnew, Xold)) /* for min heap: ( Xnew > Xold ) */
        {
            /* Xnew and Xold are out of order, so swap them */
            heap [newindex] = leaf ;
            heap [old] = i ;
        }
        else return (nheap) ;
        old = newindex ;
    }
    return (nheap) ;
}

/* ========================================================================== */
/* === napminheapify ======================================================== */
/* ========================================================================== */

/* heapify starting at node p.  On input, the heap at node p satisfies the */
/* heap property, except for heap [p] itself.  On output, the whole heap */
/* satisfies the heap property. */

PRIVATE void NAPHEAPIFY
(
    NAPINT      *heap, /* size nheap+1, containing indices into X */
    NAPFLOAT const *X, /* not modified */
    NAPINT          p, /* start at node p in the heap */
    NAPINT      nheap  /* heap [1 ... n] is in use */
)
{
    NAPFLOAT Xe, Xleft, Xright ;
    NAPINT left, right, e, hleft, hright ;

    e = heap [p] ;
    Xe = X [e] ;

    while ( TRUE )
    {
        left = p * 2 ;
        right = left + 1 ;

        if ( right <= nheap )
        {
            /* both left and right children of p are in the heap */
            hleft  = heap [left] ;
            hright = heap [right] ;
            Xleft  = X [hleft] ;
            Xright = X [hright] ;
            if (IN_ORDER (Xleft, Xright)) /* for min heap: (Xleft < Xright) */
            {
                /* traverse down the left heap */
                if (OUT_OF_ORDER (Xe, Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
            else
            {
                /* traverse down the right heap */
                if (OUT_OF_ORDER (Xe, Xright)) /* for min heap: (Xe > Xright) */
                {
                    /* swap e with the right child */
                    heap [p] = hright ;
                    p = right ;
                }
                else
                {
                    /* element e, with value Xe, has found its final home */
                    heap [p] = e ;
                    return ;
                }
            }
        }
        else
        {
            /* only the left child of p is in the heap */
            if ( left <= nheap )
            {
                /* traverse down the left heap */
                hleft = heap [left] ;
                Xleft = X [hleft] ;
                if (OUT_OF_ORDER (Xe, Xleft)) /* for min heap: (Xe > Xleft) */
                {
                    /* swap e with the left child */
                    heap [p] = hleft ;
                    p = left ;
                }
            }
            /* element e, with value Xe, has found its final home */
            heap [p] = e ;
            return ;
        }
    }
}

#undef NAPHEAP_BUILD
#undef NAPHEAP_DELETE
#undef NAPHEAP_ADD
#undef NAPHEAPIFY
#undef IN_ORDER
#undef OUT_OF_ORDER
#undef CONSTRUCT_MINHEAP_FUNCTIONS
#undef CONSTRUCT_MAXHEAP_FUNCTIONS

/* ========================================================================== */
/* === Source/napheap_check.c =============================================== */
/* ========================================================================== */

/*  Determine KKT error in a solution to the separable convex quadratic
    knapsack problem:

    min .5x'Dx - c'x  subject to  lo <= x <= hi, blo <= b <= bhi, a'x = b

    Here D is a diagonal matrix with nonnegative diagonal entries.
    Input is the solution x and the Lagrange multiplier lambda associated
    with the constraint a'x = b. The code determines the smallest perturbations
    in b, c, lo, hi, and lambda with the property that x and lambda satisfy
    the KKT conditions for the perturbed problem. The KKT conditions are the
    following:

    1. If b = a'x, then b = bhi if lambda > 0
                        b = blo if lambda < 0
                        b in [blo, bhi] if lambda = 0

    2. lo [j] <= x [j] <= hi [j]

    3. If D [j,j] = 0, then x [j] = hi [j] if a [j]*lambda - c [j] < 0
                            x [j] = lo [j] if a [j]*lambda - c [j] > 0
       If D [j,j] > 0, then

        x [j] = mid {lo [j], (c [j] - a [j]*lambda)/D [j,j], hi [j]}

       Here `mid' denotes the mean (or middle) of the three quantities.

    The code returns the relative 1-norm of the minimum perturbations
    in b, c, lo, hi, and lambda that are needed to satisfy these conditions.
    These are output in 3 arguments:

    errb: (|pert b|+|pert lambda|)/ sum |a [j]*x [j]|

    errB: ||pert lo||_1 + ||pert hi||_1 over the 1-norm of the corresponding
          original components

    errc: ||pert c||_1 over the 1-norm of the corresponding original components

    This function does not check to see if the input problem is valid or not.
    In particular, all entries in D must be >= 0, and lo [j] <= hi [j] must
    also hold for all j.  Napheap checks these conditions if Parm->check is
    true.  However, in this function, any d[j]<0 is treated as zero.
    If lo [j] > hi [j] for any j, the result of this function is undefined.
*/

NAPFLOAT napheap_check /* returns largest of the 3 error values */
(
    NAPFLOAT     *errb, /* relative error in b (only return if not NULL) */
    NAPFLOAT     *errB, /* relative error in bounds (only return if not NULL) */
    NAPFLOAT     *errc, /* relative error in c (only return if not NULL) */
    NAPdata   *napdata
)
{
    NAPINT j ;
    NAPFLOAT ax, loj, hij, t, s, err, normB, normc,
             pertb, pertB, pertc, dj ;

    NAPFLOAT const  *x = napdata->x ;      /* size n, solution */
    NAPFLOAT    lambda = napdata->lambda ; /* multiplier (scalar) */
    NAPINT           n = napdata->n ;      /* problem dimension */
    NAPFLOAT const  *c = napdata->c ;      /* size n, linear term in objective*/
    NAPFLOAT const  *d = napdata->d ;      /* size n, diagonal of cost Hessian*/
    NAPFLOAT const  *a = napdata->a ;      /* size n, linear constraint vector*/
    NAPFLOAT       blo = napdata->blo ;    /* lower bound for a'x */
    NAPFLOAT       bhi = napdata->bhi ;    /* upper bound for a'x */
    NAPparm      *Parm = napdata->Parm ;   /* parameter structure */

    int const a_is_one = (a == NULL) ? TRUE : FALSE ;

    /* get the bounds on b, and swap blo and bhi to ensure blo <= bhi */
    int const d_is_zero = (d == NULL) ? TRUE : Parm->d_is_zero ;
    int const d_is_one  = Parm->d_is_one ;
    int const c_is_zero = (c == NULL) ;
    if ( blo > bhi )
    {
        t = blo ;
        blo = bhi ;
        bhi = t ;
    }

    /* determine b perturbation */
    t = NAPZERO ;
    ax = NAPZERO ;
    for (j = 0; j < n; j++)
    {
        NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
        s = aj*x [j] ;
        t += s ;
        ax += fabs (s) ;
    }

/*printf ("blo: %e a'x: %e bhi: %e lambda: %e\n", blo, t, bhi, lambda) ;*/
    if ( lambda > NAPZERO )
    {
        pertb = fabs (t-bhi) ;
        if ( lambda < pertb ) /* less change if lambda is perturbed */
        {
            pertb = lambda ;
        }
    }
    else if ( lambda < NAPZERO )
    {
        pertb = fabs (t-blo) ;
        if ( -lambda < pertb ) /* less change if lambda is perturbed */
        {
            pertb = -lambda ;
        }
    }
    else
    {
        pertb = NAPMAX (blo-t, NAPZERO) ;
        pertb = NAPMAX (t-bhi, pertb) ;
    }
    /* error relative to sum |aj*xj| */
    if ( ax > NAPZERO ) pertb /= ax ;

    /* determine perturbations in bounds and in c */
    pertB = NAPZERO ;
    normB = NAPZERO ;
    pertc = NAPZERO ;
    normc = NAPZERO ;
    dj = NAPZERO ;
    for (j = 0; j < n; j++)
    {
        NAPFLOAT const LOj = VSZERO ;
        NAPFLOAT const HIj = VSONE ;
        loj = LOj ;
        hij = HIj ;
        if      ( x [j] < loj ) loj = x [j] ;
        else if ( x [j] > hij ) hij = x [j] ;
        dj = (NAPFLOAT) 1 ;
        NAPFLOAT const cj = c [j] ;
        if ( dj > NAPZERO )
        {
            NAPFLOAT aj = (a_is_one) ? VSONE : a [j] ;
            t = dj*x [j] + aj*lambda ; /* cj that yields x [j] */
            s = fabs (t - cj) ;           /* perturbation in c */
            /* we have the choice of perturbing c, lo, or hi */
             /*can eliminate c pert by moving hij to x [j]*/
            if ( cj > t )
            {
                if ( hij - x [j] < s ) /* pert in hi < pert in c */
                {
                    hij = x [j] ;
                    s = NAPZERO ;
                }
            }
            else /* cj <= t, can move loj to x [j] */
            {
                if ( x [j] - loj < s ) /* pert in lo < pert in c */
                {
                    loj = x [j] ;
                    s = NAPZERO ;
                }
            }
        }

        /* accumulate the errors */
        normc += fabs (cj) ;
        if ( s > NAPZERO )
        {
            pertc += s ;
        }
        if ( loj != LOj )
        {
            pertB += fabs (loj-LOj) ;
            normB += fabs (LOj) ;
        }
        if ( hij != HIj )
        {
            pertB += fabs (hij-HIj) ;
            normB += fabs (HIj) ;
        }
    }
    if ( normc > NAPZERO ) pertc /= normc ;
    if ( normB > NAPZERO ) pertB /= normB ;

    if ( errb != NULL ) *errb = pertb ;
    if ( errB != NULL ) *errB = pertB ;
    if ( errc != NULL ) *errc = pertc ;

    /* return the largest error */
    err = NAPMAX (pertb, pertc) ;
    err = NAPMAX (err, pertB) ;
    return (err) ;
}

/* Ovi's codes for finding a minimum vertex cover */
/* ========================================================================== */
// Utility Functions
/* ========================================================================== */
static VSINT *Randperm (VSINT n)
{
    VSINT *p, k, j, t ;

    p = malloc (n * sizeof (VSINT)) ;  /* allocate result */
    if (!p)
    {
        printf("Error: Not enough memory in for permutation\n") ;
        return (NULL) ;
    }
    for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;

    for (k = 0 ; k < n ; k++)
    {
        j = k + (rand ( ) % (n-k)) ;    /* j = rand integer in range k to n-1 */
        t = p [j] ;                     /* swap p[k] and p[j] */
        p [j] = p [k] ;
        p [k] = t ;
    }
    return (p) ;
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
static void Accumulate (VSLINT *p, VSLINT *c, VSINT n)
{
    VSINT i ;
    VSLINT nz = 0 ;

    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
}


static void Elapsed_Time
(
    struct timespec start, struct timespec end, VSFLOAT *elapsed
)
{
    *elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) /1e9;
}

/* ========================================================================== */
// Queue Functions
/* ========================================================================== */
// Queue is empty when ptr_end < ptr_start

static Queue *Queue_Init(VSINT size_max)
{
    Queue *Queue     = malloc(1 * sizeof(*Queue)) ;
    Queue->size_max  = size_max ;
    Queue->queue     = malloc(size_max * sizeof(*(Queue->queue))) ;
    Queue->ptr_start = 0 ;
    Queue->ptr_end   = -1 ;
    return Queue ;
}

static void Queue_Destroy(Queue *Queue)
{
    free(Queue->queue) ;
    free(Queue) ;
}

static void Queue_Push(Queue *Queue, VSINT memb)
{
    Queue->queue[ ++(Queue->ptr_end) ] = memb ;
}

static VSINT Queue_Pop(Queue *Queue)
{
    if( Queue->ptr_end < Queue->ptr_start ) return VSEMPTY ;

    VSINT memb = Queue->queue[ Queue->ptr_start ] ;
    (Queue->ptr_start)++ ;
    return memb ;
}

static CSC *CSC_Init(VSINT m, VSINT n, VSLINT nz)
{
    CSC *A = calloc(1, sizeof(*A)) ;
    A->i = calloc(nz, sizeof(*(A->i))) ;
    A->p = calloc((n+1), sizeof(*(A->p))) ;
    A->m = m ;
    A->n = n ;
    A->nz = nz ;
    return A ;
}

static void CSC_Destroy(CSC *A)
{
    free(A->i) ;
    free(A->p) ;
    free(A) ;
}

static CSC *CSC_Transpose (const CSC *A)
{
    VSINT j, *Ci, n, m, *Ai ;
    CSC *C ;
    VSLINT *Ap, *Cp, p, q, *w ;

    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ;
    C = CSC_Init(n, m, Ap[n]) ;       /* allocate result */
    w = calloc (m, sizeof (*w)) ;

    Cp = C->p ; Ci = C->i ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    Accumulate (Cp, w, m) ;                                /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
        }
    }
    free(w) ;
    return C ;
}

// Extract matrix corresponding to bipartite graph. The rows are the 0-vertices
// of the partition, and the columns are the 1-vertices
static CSC *CSC_Extract_Bipartite_Graph
(
    VSINT              *part,
    CSC                    *A,
    VSINT **localBipLabelRet
)
{
    VSINT i, j, n = A->n, *Ai = A->i, i_added, n0 = 0, n1 = 0  ;
    VSLINT end, *Ap = A->p, p, nzBip = 0 ;
    VSINT *localBipLabel = malloc(n * sizeof(*localBipLabel)) ;
    for(i = 0 ; i < n ; i++) localBipLabel[i] = VSEMPTY ;

    // Set return
    *localBipLabelRet = localBipLabel ;

    // Temporary storage for bipartite graph
    VSLINT *Bp = calloc(n+1, sizeof(*Bp)) ;
    VSINT  *Bi = malloc(Ap[n] * sizeof(*Bi)) ;

    // Build Bp and Bi.  0-vertices are columns, 1-vertices are rows
    // i runs through original graph. If it has a neighbor 
    for(i = 0 ; i < n ; i++)
    {
        if(part[i] == 0)
        {
            end = Ap[i+1] ;
            i_added = VSFALSE ;
            for(p = Ap[i] ; p < end ; p++)
            {
                j = Ai[p] ;
                if(part[j] == 1)
                {
                    if(i_added == VSFALSE)
                    {
                        // Add new vertex to 0-side of bipartite graph. Its
                        // local label is n0
                        localBipLabel[i] = n0 ;
                        i_added = VSTRUE ;
                    }

                    // Add j to the bipartite graph as a row vertex
                    if(localBipLabel[j] == VSEMPTY)
                    {
                        localBipLabel[j] = n1 ;
                        n1++ ;
                    }

                    // Add nrow to Bi as a neighbor of ncol
                    Bi[nzBip++] = localBipLabel[j] ;
                }
            }
            if(i_added)
            {
                Bp[n0 + 1] = nzBip ;
                n0++ ;
            }
        }
    }

    // At this point, the Bp and Bi arrays are correct, though Bi may not be
    // sorted. So transpose it
    CSC B ;
    B.m = n1 ; B.n = n0 ;
    B.i = Bi ; B.p = Bp ; B.nz = Bp[n0] ;
    CSC *T = CSC_Transpose(&B) ;
    free(B.i) ;
    free(B.p) ;

    return T ;
}

/* ========================================================================== */
// Min Heap Functions
/* ========================================================================== */

void MinHeap_Bubble_Down_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT   numNodes,
    VSINT  *c,
    VSINT   p
)
{
    VSINT e, Xe, hleft, hright, Xleft, Xright;
    VSINT left, right;

    e  = heap[p];
    Xe = c[e];

    while (VSTRUE)
    {
        left  = p << 1;      // p * 2
        right = left + 1;

        if (right <= numNodes)
        {
            hleft  = heap[left];
            hright = heap[right];
            Xleft  = c[hleft];
            Xright = c[hright];

            if (Xleft < Xright)
            {
                if (Xe > Xleft)
                {
                    heap[p]     = hleft;
                    lookup[hleft] = p;
                    p  = left;
                }
                else
                {
                    heap[p]   = e;
                    lookup[e] = p;
                    return;
                }
            }
            else
            {
                if (Xe > Xright)
                {
                    heap[p]      = hright;
                    lookup[hright] = p;
                    p   = right;
                }
                else
                {
                    heap[p]   = e;
                    lookup[e] = p;
                    return;
                }
            }
        }
        else
        {
            if (left <= numNodes)
            {
                hleft = heap[left];
                Xleft = c[hleft];
                if (Xe > Xleft)
                {
                    heap[p]     = hleft;
                    lookup[hleft] = p;
                    p   = left;
                }
            }
            heap[p]   = e;
            lookup[e] = p;
            return;
        }
    }
}

void MinHeap_Heapify_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT   numNodes,
    VSINT  *c
)
{
    for (VSINT p = numNodes / 2; p >= 1; p--)
    {
        MinHeap_Bubble_Down_v0(heap, lookup, numNodes, c, p);
    }
}

void MinHeap_Insert_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT   toInsert
)
{
    VSINT node_child = ++(*numNodes);
    VSINT c_child    = c[toInsert];
    VSINT node_par;

    while (node_child > 1)
    {
        node_par = node_child >> 1; // node_child / 2
        if (c[heap[node_par]] > c_child)
        {
            heap[node_child]           = heap[node_par];
            lookup[heap[node_par]]     = node_child;
            node_child                 = node_par;
        }
        else
        {
            break;
        }
    }

    heap[node_child]   = toInsert;
    lookup[toInsert]   = node_child;
}

// ASSUMED: toDel is in the heap
void MinHeap_Delete_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT   toDel
)
{
    VSINT address = lookup[toDel];
    if (address == VSEMPTY) return;
    if (address == *numNodes)
    {
        (*numNodes)--;
        return;
    }

    VSINT replacement = heap[(*numNodes)--];
    lookup[toDel]      = VSEMPTY;

    if (address == 1 || c[heap[address >> 1]] <= c[replacement])
    {
        heap[address]    = replacement;
        lookup[replacement] = address;
        MinHeap_Bubble_Down_v0(heap, lookup, *numNodes, c, address);
    }
    else
    {
        while (address > 1)
        {
            VSINT par = address >> 1;
            if (c[heap[par]] > c[replacement])
            {
                heap[address]           = heap[par];
                lookup[heap[par]]       = address;
                address                 = par;
            }
            else
            {
                break;
            }
        }
        heap[address]     = replacement;
        lookup[replacement] = address;
    }
}

void MinHeap_See_Top_Val_v0
(
    VSINT  *heap,
    VSINT  *c,
    VSINT   numNodes,
    VSINT  *topVal
)
{
    *topVal = (numNodes != 0) ? heap[1] : VSEMPTY;
}

void MinHeap_Get_Top_Val_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT  *topVal
)
{
    if (*numNodes == 0)
    {
        *topVal = VSEMPTY;
        return;
    }
    *topVal = heap[1];
    MinHeap_Delete_v0(heap, lookup, numNodes, c, *topVal);
}

void MinHeap_See_Top_Key_v0
(
    VSINT  *heap,
    VSINT  *c,
    VSINT   numNodes,
    VSINT  *topKey
)
{
    *topKey = (numNodes != 0) ? c[heap[1]] : VSEMPTY;
}

void MinHeap_Get_Top_Key_v0
(
    VSINT  *heap,
    VSINT  *lookup,
    VSINT  *numNodes,
    VSINT  *c,
    VSINT  *topKey
)
{
    if (*numNodes == 0)
    {
        *topKey = VSEMPTY;
        return;
    }
    *topKey = c[heap[1]];
    MinHeap_Delete_v0(heap, lookup, numNodes, c, heap[1]);
}


/* ========================================================================== */
// Max Matching Function
/* ========================================================================== */
/* find an augmenting path starting at column k and extend the match if found */
static void Augment (VSINT k, const CSC *A, VSINT *jmatch, VSLINT *cheap,
        VSINT *w, VSINT *js, VSINT *is, VSLINT *ps)
{
    VSLINT *Ap = A->p, p ;
    VSINT found = 0, i = -1, *Ai = A->i, head = 0, j ;
    js [0] = k ;                        /* start with just node k in jstack */
    while (head >= 0)
    {
        /* --- Start (or continue) depth-first-search at node j ------------- */
        j = js [head] ;                 /* get j from top of jstack */
        if (w [j] != k)                 /* 1st time j visited for kth path */
        {
            w [j] = k ;                 /* mark j as visited for kth path */
            for (p = cheap [j] ; p < Ap [j+1] && !found ; p++)
            {
                i = Ai [p] ;            /* try a cheap assignment (i,j) */
                found = (jmatch [i] == -1) ;
            }
            cheap [j] = p ;             /* start here next time j is traversed*/
            if (found)
            {
                is [head] = i ;         /* column j matched with row i */
                break ;                 /* end of augmenting path */
            }
            ps [head] = Ap [j] ;        /* no cheap match: start dfs for j */
        }
        /* --- Depth-first-search of neighbors of j ------------------------- */
        for (p = ps [head] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;                /* consider row i */
            if (w [jmatch [i]] == k) continue ; /* skip jmatch [i] if marked */
            ps [head] = p + 1 ;         /* pause dfs of node j */
            is [head] = i ;             /* i will be matched with j if found */
            js [++head] = jmatch [i] ;  /* start dfs at column jmatch [i] */
            break ;
        }
        if (p == Ap [j+1]) head-- ;     /* node j is done; pop from stack */
    }                                   /* augment the match if path found: */
    if (found) for (p = head ; p >= 0 ; p--) jmatch [is [p]] = js [p] ;
}

/* find a maximum transversal */
 /*[jmatch [0..m-1]; imatch [0..n-1]]*/
static VSINT *Maximum_Matching_Get(const CSC *A, VSINT seed) 
{
    VSINT i, j, k, n, m, n2 = 0, m2 = 0, *jimatch,
    *Ai, *jmatch, *imatch, *q, *w ;
    CSC *C ;
    VSLINT *Ap, *Cp, p ;

    n = A->n ; m = A->m ; Ap = A->p ; Ai = A->i ;
    w = jimatch = calloc (m+n, sizeof (VSINT)) ;   /* allocate result */

    for (k = 0, j = 0 ; j < n ; j++)    /* count nonempty rows and columns */
    {
        n2 += (Ap [j] < Ap [j+1]) ;
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            w [Ai [p]] = 1 ;
            k += (j == Ai [p]) ;        /* count entries already on diagonal */
        }
    }
    if (k == VSMIN (m,n))              /* quick return if diagonal zero-free */
    {
        jmatch = jimatch ; imatch = jimatch + m ;
        for (i = 0 ; i < k ; i++) jmatch [i] = i ;
        for (      ; i < m ; i++) jmatch [i] = -1 ;
        for (j = 0 ; j < k ; j++) imatch [j] = j ;
        for (      ; j < n ; j++) imatch [j] = -1 ;
        return jimatch ;
    }
    for (i = 0 ; i < m ; i++) m2 += w [i] ;
    C = (m2 < n2) ? CSC_Transpose (A) : ((CSC *) A) ; /* transpose if needed */

    n = C->n ; m = C->m ; Cp = C->p ;
    jmatch = (m2 < n2) ? jimatch + n : jimatch ;
    imatch = (m2 < n2) ? jimatch : jimatch + m ;

    // w = malloc(5*n * sizeof (VSINT)) ;             /* get workspace */

    // cheap = w + n ; js = w + 2*n ; is = w + 3*n ; ps = w + 4*n ;

    w  = malloc(n * sizeof(*w)) ;
    VSLINT *cheap = malloc(n * sizeof(*cheap)) ;
    VSINT  *js = malloc(n * sizeof(*js)) ;
    VSINT  *is = malloc(n * sizeof(*is)) ;
    VSLINT *ps = malloc(n * sizeof(*ps)) ;


    for (j = 0 ; j < n ; j++) cheap [j] = Cp [j] ; /* for cheap assignment */
    for (j = 0 ; j < n ; j++) w [j] = -1 ;         /* all columns unflagged */
    for (i = 0 ; i < m ; i++) jmatch [i] = -1 ;    /* nothing matched yet */
    q = Randperm (n) ;                   /* q = random permutation */
    for (k = 0 ; k < n ; k++)   /* augment, starting at column q[k] */
    {
        Augment (q ? q [k]: k, C, jmatch, cheap, w, js, is, ps) ;
    }
    free(q) ;
    for (j = 0 ; j < n ; j++) imatch [j] = -1 ;     /* find row match */
    for (i = 0 ; i < m ; i++) if (jmatch [i] >= 0) imatch [jmatch [i]] = i ;
    if( (m2 < n2) ) CSC_Destroy(C) ;
    free(w) ;
    free(cheap) ;
    free(js) ;
    free(is) ;
    free(ps) ;
    return jimatch ;
}

/* ========================================================================== */
// Min Vertex Cover Function
/* ========================================================================== */
// The graph is given in CSC format
// The two sides of the bipartite graphs are cols and rows. Each set has a local
// indexing, cols are 0...n-1 and rows are 0...m-1
// The global indexing is rows then cols. So local indexing coincides with
// global indexing for rows, and global index of a row is its local index
// shiften by m
// Matching: at first it's a list of the matches of row indices, then the
// matches of col indices. In particular, matching[0] is the match of row index
// 0 (or VSEMPTY if row index 0 is not matched)
static void Vtx_Cover_From_Max_Matching_Get
(
    CSC            *B,
    VSINT   *matching,
    VSINT     **cover,
    int    start_from_rows,   // NEW
    CSC            *Bt        // NEW: required iff start_from_rows
)
{
    VSINT m = B->m, n = B->n, i, j, k;
    VSLINT p, end, *Sp;
    VSINT  *Si;
    VSINT  *matchingStart, *matchingOther;
    VSINT  *coverStart,    *coverOther;
    VSINT  n_start, n_other;

    if (!start_from_rows) {
        Sp = B->p;  Si = B->i;
        n_start = n;  n_other = m;
        matchingStart = matching + m;   // col matches
        matchingOther = matching;        // row matches
    } else {
        Sp = Bt->p; Si = Bt->i;
        n_start = m;  n_other = n;
        matchingStart = matching;        // row matches
        matchingOther = matching + m;    // col matches
    }

    *cover = malloc((m + n) * sizeof(**cover));
    if (!start_from_rows) {
        coverStart = (*cover) + m;       // col cover slot
        coverOther = *cover;              // row cover slot
    } else {
        coverStart = *cover;              // row cover slot
        coverOther = (*cover) + m;       // col cover slot
    }

    for (i = 0; i < n_other; i++) coverOther[i] = 0;
    for (i = 0; i < n_start; i++) coverStart[i] = 1;

    VSINT *visited = calloc(m + n, sizeof(*visited));
    Queue *Q = Queue_Init(m + n);

    // Seed the queue with unmatched vertices on the starting side.
    // Queue convention: [0, n_other) are "other" verts,
    //                   [n_other, n_other + n_start) are "starting" verts.
    for (i = 0; i < n_start; i++) {
        if (matchingStart[i] != VSEMPTY) continue;
        coverStart[i] = 0;
        visited[n_other + i] = 1;
        Queue_Push(Q, n_other + i);
    }

    while ((j = Queue_Pop(Q)) != VSEMPTY) {
        if (j < n_other) {
            coverOther[j] = 1;
            k = matchingOther[j];
            if (k != VSEMPTY && !visited[n_other + k]) {
                Queue_Push(Q, n_other + k);
                visited[n_other + k] = 1;
            }
        } else {
            j = j - n_other;
            coverStart[j] = 0;
            end = Sp[j+1];
            for (p = Sp[j]; p < end; p++) {
                i = Si[p];
                if (visited[i]) continue;
                if (matchingStart[j] != i) {
                    Queue_Push(Q, i);
                    visited[i] = 1;
                }
            }
        }
    }

    Queue_Destroy(Q);
    free(visited);
}

// Check if all neighbors of vtx_col are in the cover. It is 
static int Are_All_Neighbors_In_Cover
(
    CSC             *B,
    VSINT     vtx_col,
    VSINT      *cover  // cover[i] = 1 if row vertex i is in the cover  
)
{
    VSLINT end, p, *Bp = B->p ;
    VSINT i, *Bi = B->i ;

    end = Bp[vtx_col+1] ;
    for(p = Bp[vtx_col] ; p < end ; p++)
    {
        i = Bi[p] ;
        if(cover[i] == 0)
        {
            return VSFALSE ;
            break ;
        }
    }
    return VSTRUE ;
}

// Thm: a cover is NOT minimal iff there exists a node in the cover s.t. all its
// neighbors are in the cover.
// Proof: Suppose the cover is not minimal. Then there is a node you can remove
// and still retain a cover. This node you can remove, all its neighbors must be
// in the cover, otherwise its removal would break the cover. So done.
// Suppose there is a node in the cover s.t. all its neighbors are in the cover
// too. Then you can remove it and still retain a cover, so it was not minimal.
// Done.
static void Check_Bip_Vertex_Cover(CSC *B, VSINT *cover)
{
    VSINT *coverRows = cover ;
    VSINT *coverCols = cover + B->m ;

    VSINT m = B->m, n = B->n, j ;

    // Check that at least one endpoint of each edge is in the cover
    // If j is not in the cover, all its neighbors i HAVE to be in the cover,
    // otherwise edge (i, j) is not covered. This also serves to show
    // minimality. If a col vertex is not in the cover, then all its neighbors
    // must be in the cover, and none of them can be removed. If a
    for(j = 0 ; j < n ; j++)
    {
        if(cover[m+j] == 0)
        {
            if(!Are_All_Neighbors_In_Cover(B, j, coverRows))
            {
                printf
                (
                    "Error: Col vertex %lld and one of its neighbors are not"
                    "in the cover, meaning it's not a cover", PRI(j)
                ) ;
                break ;
            }
        }
        else // j is in the cover
        {
            // Not all of js neighbors should be in the cover. If they are,
            // the cover is not minimal as j is redundant
            if(Are_All_Neighbors_In_Cover(B, j, coverRows))
            {
                printf
                (
                    "Error: All neighbors of col vtx %lld are in the cover"
                    ", so the cover is not minimal\n", PRI(j)
                ) ;
                break ;
            }
        }
    }

    // Transpose
    CSC *T = CSC_Transpose(B) ;

    // Do the same for row vertices of the original matrix (col vertices of T)
    for(j = 0 ; j < T->n ; j++)
    {
        if(cover[j] == 0)
        {
            if(!Are_All_Neighbors_In_Cover(T, j, coverCols))
            {
                printf
                (
                    "Error: Row vertex %lld and one of its neighbors are not"
                    "in the cover, meaning it's not a cover", PRI(j)
                ) ;
                break ;
            }
        }
        else // j is in the cover
        {
            // Not all of js neighbors should be in the cover. If they are,
            // the cover is not minimal as j is redundant
            if(Are_All_Neighbors_In_Cover(T, j, coverCols))
            {
                printf
                (
                    "Error: All neighbors of row vtx %lld are in the cover,"
                    ", so the cover is not minimal\n", PRI(j)
                ) ;
                break ;
            }
        }
    }
    CSC_Destroy(T) ;
}

static void Check_Matching_Size_Against_Cover_Size
(
    CSC *B,
    VSINT *matching,
    VSINT *cover
)
{
    // Get matching size
    VSINT matchingSize = 0, numRows = B->m, coverSize = 0,
    numNodes = B->m + B->n, i ;

    for(i = 0 ; i < numRows ; i++)
    {
        if(matching[i] != VSEMPTY) matchingSize++ ;
    }

    // Get cover size
    for(i = 0 ; i < numNodes ; i++)
    {
        if(cover[i]) coverSize++ ;
    }

    if(coverSize != matchingSize)
    {
        printf
        (
            "Error: cover size = %lld doesn't match matching size = %lld\n",
            PRI(coverSize), PRI(matchingSize)
        ) ;
    }
}

/* ========================================================================== */
// Separator Functions
/* ========================================================================== */

// Print error message if we do NOT have a separator
// Silent if the partition is a separator
static void Separator_Check
(
    VSINT *part, // entries are 0, 1, 2. Part 2 is candidate separator
    CSC       *A  // Adjacency matrix
)
{
    VSINT i, j, *Ai = A->i, n = A->n ;
    VSLINT p, *Ap = A->p, end ;
    for(i = 0 ; i < n ; i++)
    {
        if(part[i] == 0)
        {
            end = Ap[i+1] ;
            for(p = Ap[i] ; p < end ; p++)
            {
                j = Ai[p] ;
                if(part[j] == 1)
                {
                    printf
                    (
                        "\nError: vertex %lld in partition 0 and vertex %lld "
                        "in partition 1\n", (long long) i, (long long) j   
                    ) ;
                    return ;
                }
            }
        }
    }
}

// Given a graph and a partition part, make part 2 a separator by finding a
// vertex cover of the bipartite graph corresponding to edges cut by parts
// 0, 1 and moving the vertices of the cover into part 2. Modify the part array
// in-place.
// Inputs:
// - part: array of length nwith entries 0, 1, or 2 signifying which partition
// a vertex is in
// - n, Ai, Ap: number of vertices, and adjacency lists and column pointers
// storing the graph
// - vtx_wgts: array of length n storing vertex weights
// - get_stats: a flag which when enabled, signals the function to computte and
// return various stats
// - n_bip, nz_bip: holds the number of vertices in the bipartite graph
// and the number of non-zeros in its matrix representation
// - wgt_cover: holds the weight of the cover
// - imbal: abs. value of difference in vertex weight between the two sides
// of the bipartite graph, ignoring the vertices in the vertex cover
// - time_gen_bip_graph_ret: holds the time required to generate the matrix
// representation of the bipartite graph. Note: for this algorithm, vertex
// weights of the bipartite graph are not used, hence they don't need to be
// computed
// - time_cover_ret: holds the time required just to compute a vertex cover
// of the bipartite graph


// Inner workings notes:
// We distinguish between 0-vertices and 1-vertices
// localBipLabel[i]: VSEMPTY if i is not in the bipartite graph. Otherwise, if
//                   part[i] == 0 then it's the local index in the 0-side,
//                   or if part[i] == 1 then it's the local index in the 1-side
// 
// As far as the global label in the bipartite graph, row vertices are labeled
// first, and then column vertices are labeled second 
void vs_separator
(
    int                *part, // array with entries 0, 1, or 2
    VSINT                *Ai, // Original graph
    VSLINT               *Ap, // Original graph
    VSINT                  n, // Original graph n
    VSINT          *vtx_wgts,
    int      try_both_covers, // If nonzero: compute both Konig covers (BFS
                              // from unmatched cols and from unmatched rows)
                              // and pick the one with smallest weight,
                              // tiebreaking on smallest imbalance.
                              // If zero: original logic, single cover from
                              // BFS on unmatched columns.
    int            get_stats,  // Binary flag to determine if (somewhat
                               // expensive) stats should be computed

    // Ignored if get_stats is false. In particular, they can be NULL
    VSINT             *n_bip,
    VSLINT           *nz_bip,
    VSINT         *wgt_cover,
    VSINT             *imbal,
    VSFLOAT *time_entire_ret
)
{
    // The original (single-cover) stats path computes its stats from a
    // bipartite weight array w[] AFTER part[] has been overwritten with
    // separators, so it needs a copy of the original partition. The
    // two-cover path doesn't (it computes everything before modifying part[]).
    VSINT *part_org_cpy = NULL ;
    if(get_stats && !try_both_covers)
    {
        part_org_cpy = malloc(n * sizeof(*part_org_cpy)) ;
        for(VSINT i = 0 ; i < n ; i++)
        {
            part_org_cpy[i] = part[i] ;
        }
    }
    struct timespec start, end, start_entire, end_entire ;
    if(get_stats) TIME(start_entire) ;
    CSC A ;
    A.m = n, A.n = n ; A.nz = Ap[n] ; A.i = Ai ; A.p = Ap ;
    VSINT *localBipLabel, i, j, *matching ;
    VSFLOAT time_generate, time_cover, time_cover_to_sep, time_entire ;

    // This alg doesn't use the vertex weights of the bipartite graph
    // So they do not need to be genererated, and hence do not need
    // to be timed
    TIME(start) ;
    CSC *B = CSC_Extract_Bipartite_Graph(part, &A, &localBipLabel) ;
    TIME(end) ;
    Elapsed_Time(start, end, &time_generate) ;

    VSINT nrows = B->m ;

    // Get maximum matching
    TIME(start) ;
    matching = Maximum_Matching_Get(B, 1) ;

    VSINT *cover ;
    VSINT chosen_wgt = 0, chosen_imbal = 0 ;  // only used in try_both_covers

    if(try_both_covers)
    {
        // Two candidate Konig covers (same cardinality, possibly different
        // weight and balance). Heuristic for min-weight vertex cover: pick
        // the cheaper of the two; break ties on partition balance.
        //   cover_col_start: BFS from unmatched columns (standard)
        //   cover_row_start: BFS from unmatched rows (needs B^T)
        CSC *Bt = CSC_Transpose(B) ;
        VSINT *cover_col_start, *cover_row_start ;
        Vtx_Cover_From_Max_Matching_Get(B, matching, &cover_col_start, 0, NULL);
        Vtx_Cover_From_Max_Matching_Get(B, matching, &cover_row_start, 1, Bt)  ;
        CSC_Destroy(Bt) ;

        // For each cover: cover weight (primary, minimize) and imbalance
        // (tiebreak, minimize). part[] still holds original 0/1 values here.
        VSINT wgt_cov_a = 0, wgt0_a = 0, wgt1_a = 0 ;
        VSINT wgt_cov_b = 0, wgt0_b = 0, wgt1_b = 0 ;
        for(i = 0 ; i < n ; i++)
        {
            j = localBipLabel[i] ;
            if(j == VSEMPTY) continue ;
            VSINT g  = (part[i] == 0) ? j : nrows + j ;
            VSINT wi = vtx_wgts[i] ;

            if(cover_col_start[g])      wgt_cov_a += wi ;
            else if(part[i] == 0)       wgt0_a    += wi ;
            else                        wgt1_a    += wi ;

            if(cover_row_start[g])      wgt_cov_b += wi ;
            else if(part[i] == 0)       wgt0_b    += wi ;
            else                        wgt1_b    += wi ;
        }
        VSINT imbal_a = VSABS(wgt0_a - wgt1_a) ;
        VSINT imbal_b = VSABS(wgt0_b - wgt1_b) ;

        int pick_a ;
        if      (wgt_cov_a < wgt_cov_b) pick_a = 1 ;
        else if (wgt_cov_a > wgt_cov_b) pick_a = 0 ;
        else                            pick_a = (imbal_a <= imbal_b) ;

        if(pick_a)
        {
            cover        = cover_col_start ;
            chosen_wgt   = wgt_cov_a ;
            chosen_imbal = imbal_a ;
            free(cover_row_start) ;
        }
        else
        {
            cover        = cover_row_start ;
            chosen_wgt   = wgt_cov_b ;
            chosen_imbal = imbal_b ;
            free(cover_col_start) ;
        }
        // printf("%s,%i,%i,%i,%i\n", G_Current_Problem_Name, wgt_cov_a, wgt_cov_b, imbal_a, imbal_b) ;
        fflush(stdout) ;
    }
    else
    {
        // Original logic: single Konig cover from BFS on unmatched columns.
        Vtx_Cover_From_Max_Matching_Get(B, matching, &cover, 0, NULL) ;
    }
    TIME(end) ;
    Elapsed_Time(start, end, &time_cover) ;


    // Move vertices in the cover to separator
    TIME(start) ;
    VSINT offset ;
    for(i = 0 ; i < n ; i++)
    {
        j = localBipLabel[i] ;
        if(j != VSEMPTY)
        {
            offset = ( part[i] ) ? nrows : 0 ;
            if(cover[offset + j]) part[i] = 2 ;
        }
    }
    TIME(end) ;
    Elapsed_Time(start, end, &time_cover_to_sep) ;
    // printf("\nCover to sep: %f\n", time_cover_to_sep) ;



    // Compute stats
    if(get_stats)
    {
        if(try_both_covers)
        {
            // Values already computed during cover selection
            TIME(end_entire) ;
            Elapsed_Time(start_entire, end_entire, time_entire_ret) ;

            *n_bip     = B->m + B->n ;
            *nz_bip    = B->nz ;
            *wgt_cover = chosen_wgt ;
            *imbal     = chosen_imbal ;
        }
        else
        {
            // Original stats path: build bipartite vertex weights, then
            // compute cover weight and imbalance from them.
            VSINT *w = malloc( (B->m + B->n) * sizeof(*w) ) ;
            for(i = 0 ; i < n ; i++)
            {
                j = localBipLabel[i] ;
                if(j == VSEMPTY) continue ;

                if(part_org_cpy[i] == 0) // i is a row vertex
                {
                    w[j] = vtx_wgts[i] ;
                }
                else if(part_org_cpy[i] == 1) // i is a column vertex
                {
                    w[nrows+j] = vtx_wgts[i] ;
                }
            }

            // Get the time
            TIME(end_entire) ;
            Elapsed_Time(start_entire, end_entire, time_entire_ret) ;

            TIME(start) ;

            *n_bip  = B->m + B->n ;
            *nz_bip = B->nz ;

            // Return weight of cover. Loop through all vertices of the graph
            *wgt_cover = 0 ;
            for(i = 0 ; i < B->m + B->n ; i++)
            {
                *wgt_cover += cover[i] * w[i] ;
            }


            // Return imbalance of cover
            // Compute weight of rows vertices not in cover
            VSINT wgt_row = 0, wgt_col = 0 ;
            for(i = 0 ; i < B->m + B->n ; i++)
            {
                if(cover[i]) continue ;

                if(i < B->m) wgt_row += w[i] ;
                else         wgt_col += w[i] ;
            }

            *imbal = VSABS(wgt_row - wgt_col) ;
            free(w) ;
            free(part_org_cpy) ;

            TIME(end) ;
            VSFLOAT time_stats ;
            Elapsed_Time(start, end, &time_stats) ;
            // printf("\nStats time: %f\n", time_stats) ;
        }
    }

    #ifdef DEBUG
        Check_Bip_Vertex_Cover(B, cover) ;
        Check_Matching_Size_Against_Cover_Size(B, matching, cover) ;
        Separator_Check(part, &A) ;
    #endif

    free(cover) ;
    free(matching) ;
    free(localBipLabel) ;

    CSC_Destroy(B) ;
}

static void Check_Its_A_Cover(CSC *B, VSINT *cover)
{
    VSINT i, j, n = B->n, m = B->m, *Bi = B->i ;
    VSLINT p, *Bp = B->p ;
    for(j = 0 ; j < n ; j++)
    {
        for(p = Bp[j] ; p < Bp[j+1] ; p++)
        {
            i = Bi[p] ;
            if(!(cover[i]) && !(cover[m+j]))
            {
                printf
                (
                    "Error: Edge {%lld, %lld} is not covered by the vertex "
                    "cover\n", PRI(i), PRI(j)
                ) ;
                goto DONE ;
            } 
        }
    }
    DONE: ;
}

void Pricing_Algorithm_Vertex_Cover_Get
(
    CSC            *B,
    VSINT         *w,
    VSINT     *cover
)
{
    VSLINT p, end, *Bp = B->p ;
    VSINT i, j, *Bi = B->i, n = B->n, m = B->m, *budget, j_global, price ;

    budget = malloc((m + n) * sizeof(*budget)) ;
    for(i = 0 ; i < m + n ; i++) budget[i] = w[i] ;

    // Let i be a global vertex index. Then budget[i] = weight of i
    // which is currently used 
    for(j = 0 ; j < n ; j++) // j is a local column index
    {
        if(budget[m + j] == 0) continue ;

        j_global = m + j ;
        end = Bp[j+1] ;

        p = Bp[j] ;
        while(p < end && budget[j_global])
        {
            i = Bi[p++] ;
            price = VSMIN(budget[i], budget[j_global]) ;
            budget[i]        -= price ;
            budget[j_global] -= price ;
        }
    }

    // Set the cover according to the budget
    VSINT num_vert = m + n ;
    for(i = 0 ; i < num_vert ; i++)
    {
        cover[i] = !(budget[i]) ;
    }

    free(budget) ;
}

// If get_stats is true, return some statistics about the bipartite graph
// and vertex cover
void vs_separator_approx
(
    int                *part, // array with entries 0, 1, or 2
    VSINT                *Ai, // Original graph
    VSLINT               *Ap, // Original graph
    VSINT                  n, // Original graph n
    VSINT          *vtx_wgts,
    int            get_stats,  // Binary flag to determine if (somewhat
                               // expensive) stats should be computed

    // Ignored if get_stats is false. In particular, they can be NULL
    VSINT       *num_vert_bip,
    VSLINT            *nz_bip,
    VSINT          *wgt_cover,
    VSINT              *imbal,
    VSFLOAT *time_entire_ret
)
{
    struct timespec start, end, start_entire, end_entire ;
    if(get_stats) TIME(start_entire) ;
    CSC A ;
    VSFLOAT elapsed_generate, elapsed_cover ;
    A.m = n, A.n = n ; A.nz = Ap[n] ; A.i = Ai ; A.p = Ap ;
    VSINT *localBipLabel, i, j, m, offset ;

    // TODO: I think that paradoxically, generating w takes quite
    // a lot of time. 
    clock_gettime(CLOCK_MONOTONIC, &start) ;
    CSC *B = CSC_Extract_Bipartite_Graph(part, &A, &localBipLabel) ;

    VSINT ncols = B->n, nrows = B->m ;
    VSINT *w = malloc((nrows + ncols) * sizeof(*w)) ;
    for(i = 0 ; i < n ; i++)
    {
        j = localBipLabel[i] ;
        if(j == VSEMPTY) continue ;

        if(part[i] == 0) // i is a row vertex
        {
            w[j] = vtx_wgts[i] ;
        }
        else if(part[i] == 1) // i is a column vertex
        {
            w[nrows+j] = vtx_wgts[i] ; 
        }
    }
    
    VSINT *cover = malloc((B->m + B->n) * sizeof(*cover) ) ;

    clock_gettime(CLOCK_MONOTONIC, &end) ;
    Elapsed_Time(start, end, &elapsed_generate) ;

    // Time getting the cover
    clock_gettime(CLOCK_MONOTONIC, &start);
    Pricing_Algorithm_Vertex_Cover_Get(B, w, cover) ;
    clock_gettime(CLOCK_MONOTONIC, &end);
    Elapsed_Time(start, end, &elapsed_cover) ;

    // Compute stats
    if(get_stats)   
    {
        TIME(end_entire) ;
        Elapsed_Time(start_entire, end_entire, time_entire_ret) ; 

        *num_vert_bip = B->m + B->n ;
        *nz_bip = B->nz ;

        // Return weight of cover. Loop through all vertices of the graph
        *wgt_cover = 0 ;
        for(i = 0 ; i < B->m + B->n ; i++)
        {
            *wgt_cover += cover[i] * w[i] ;
        }


        // Return imbalance of cover
        // Compute weight of rows vertices not in cover
        VSINT wgt_row = 0, wgt_col = 0 ;
        for(i = 0 ; i < B->m + B->n ; i++)
        {
            if(cover[i]) continue ;

            if(i < B->m) wgt_row += w[i] ;
            else         wgt_col += w[i] ;
        }

        *imbal = VSABS(wgt_row - wgt_col) ;
    }

    // Move vertices in the cover to separator
    // TODO: This seems to take a long time too
    for(i = 0 ; i < n ; i++)
    {
        j = localBipLabel[i] ;
        if(j != VSEMPTY)
        {
            offset = ( part[i] ) ? nrows : 0 ;
            if(cover[offset + j]) part[i] = 2 ;
        }
    }

    #ifdef DEBUG
        // Check_Bip_Vertex_Cover(B, cover) ;
        Check_Its_A_Cover(B,cover) ;
        Separator_Check(part, &A) ;
    #endif

    free(cover) ;
    free(localBipLabel) ;
    free(w) ;
    CSC_Destroy(B) ;
}

// Update gy, the gradient with respect to y, as a result of x_i changing by
// delta_xi:
// new_gy = old_gy - (delta_xi) * gamma * A_i,
// where A_i is column of A corresponding to index i
static void Update_Gradient
(
    VSFLOAT   delta_xi,
    VSINT     i,
    VSFLOAT  *gy,
    VSINT    *Ai,
    VSLINT   *Ap,
    VSFLOAT   gamma
)
{
    VSLINT p, end = Ap[i+1] ;
    VSFLOAT const coeff = delta_xi * gamma ;
    for(p = Ap[i] ; p < end ; p++)
    {
        gy[ Ai[p] ] -= coeff ;
    }
    gy[i] -= coeff ;
}

static VSINT VSP_interior_and_penalty_reduction
(
   VSINT    ilower, /* search x for indices > ilower and < n */
   VSFLOAT      *x, /* search for 0 < x [i] < 1 */
   VSFLOAT     *gx,
   VSFLOAT     *gy,
   VSINT         n, /* dimension of x */
   VSINT       *Ai,
   VSLINT      *Ap,
   VSFLOAT   gamma
)
{
    VSINT i ;
    VSFLOAT xi ;

    for (i = ilower+1; i < n; i++)
    {
        xi = x[i] ;

        // Set x_i to 0 if doing so doesn't damage the objective
        if(xi > 0 && gx[i] <= VSZERO_PLUS)
        {
            x[i] = 0 ;
            Update_Gradient(-xi, i, gy, Ai, Ap, gamma) ;
        }

        // 
        if ( (xi > VSZERO) && (xi < VSONE) ) break ;
    }
    if ( i < n ) return (i) ;
    else         return (VSEMPTY) ;
}

/* =========================================================================
   ============================= VSPboundary ===============================
   =========================================================================
       ________________________________________________________________
      |    Consider the following problem:                             |
      |                                                                |
      |    (P)        max f(x) := x'g                                  |
      |                                                                |
      |               subject to a'x = m, 0 <= x <= 1                  |
      |                                                                |
      |    n = number of vertices in an undirected graph, a > 0, m     |
      |    arbitrary integer in [0, a'*1]. Given a feasible point x,   |
      |    the code replaces x with a new feasible point x* satisfying |
      |    f(x*) >= f(x).  Moreover, x [i] = 0 or 1 for at least n-1   |
      |    distinct i in [0, n). If x* is binary, then the code returns|
      |    -1. If x* is fractional, then the code returns the index of |
      |    the fractional components of x*.                            | 
      |                                                                |
      |    Algorithm:                                                  |
      |                                                                |
      |    First search for i and then search for j such that          |
      |                                                                |
      |        0 < x [i] < 1, 0 < x [j] < 1, and i < j.                |
      |                                                                |
      |    If i does not exist, then return -1. If j does not exist,   |
      |    then return j. Assuming both i and j exist, then let        |
      |    p denote a vector with all components zero except for       |
      |    p [i] = a [j] and p [j] = -a [i]. For this choice of p,     |
      |    we have a'p = 0. Hence, x(s) := x + s*p satisfies the linear|
      |    constraint of (P) for any choice of s. By the structure of  |
      |    the objective, we have                                      |
      |                                                                |
      |        f(x + s*p) = f(x) + s*b  where  b := g_i p_i + g_j p_j. |
      |                                                                |
      |    If b >= 0, then let s > 0 be the first value of s > 0 where |
      |    either x_i (s) = 1  or x_j (s) = 0. If b < 0, then let s < 0|
      |    be the first value of s < 0 where either x_i (s) = 0 or     |
      |    x_j (x) = 1. If s is choosen in this way, the objective     | 
      |    value cannot decrease and the new x remains feasible.       |
      |________________________________________________________________|

       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|*/

VSINT VSPboundary
(
    VSFLOAT           *x, /* feasible point for (GB), output described above */
    VSFLOAT           *a, /* vertex weights, NULL => all entries in a are 1 */
    VSFLOAT          *gx, /* linear vector in objective */
    VSFLOAT          *gy,
    VSINT              n, /* problem dimension */
    VSINT            *Ai,
    VSLINT           *Ap,
    VSFLOAT        gamma,
    int const PrintLevel
)
{
    int dumpi, dumpj, two_hit ;
    VSFLOAT pi, pj, s1, s2, s3, s4, xi, xj, tlinear ;
    VSINT i, j ;

    /* find index > -1 */
    i = VSP_interior_and_penalty_reduction(-1, x, gx, gy, n, Ai, Ap, gamma) ;
    if ( i == VSEMPTY )
    {
        /* x has no interior components */
        return (EMPTY) ;
    }
    j = VSP_interior_and_penalty_reduction(i, x, gx, gy, n, Ai, Ap, gamma) ;
    if ( j == VSEMPTY )
    {
        /* done since x has at most one fractional component */
        return (i) ;
    }
    int const a_is_one = (a == NULL) ? VSTRUE : VSFALSE ;
    two_hit = VSFALSE ;
    dumpi = VSFALSE ;
    dumpj = VSFALSE ;
    while ( 1 )
    {
        /* components i and j of x are fractional */
        VSFLOAT const ai = (a_is_one) ? VSONE : a [i] ;
        VSFLOAT const aj = (a_is_one) ? VSONE : a [j] ;
        VSFLOAT const xi = x [i] ;
        VSFLOAT const xj = x [j] ;
        /* pi = aj, pj = -ai */
        tlinear = gx [i]*aj - gx [j]*ai ;

        /* if tlinear = 0, then determine whether a positive or a negative s
           step is larger in absolute value, and take the larger step */
        if ( tlinear >= VSZERO )
        {
            /* pos s step */
            s1 = (1-xi)/aj ; /* xi hits upper bound */
            s2 = xj/ai ;     /* xj hits lower bound */
            if ( s1 < s2 ) /* s is positive and xi hits upper bound */
            {
                pi = VSONE - xi ;
                x [i] = VSONE ;
                Update_Gradient(pi, i, gy, Ai, Ap, gamma) ;
                pj = -s1*ai ;
                x [j] += pj ;
                Update_Gradient(pj, j, gy, Ai, Ap, gamma) ;
                if ( x [j] <= VSZERO )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %25.15e violates lower bound\n",
                                 j, x [j]) ;
                    }
                    two_hit = VSTRUE ;
                    x [j] = VSZERO ;
                }
                else dumpi = VSTRUE ;
            }
            else /* xj hits lower bound */
            {
                pj = -xj ;
                x [j] = VSZERO ;
                Update_Gradient(pj, j, gy, Ai, Ap, gamma) ;
                pi = s2*aj ;
                x [i] += pi ;
                Update_Gradient(pi, i, gy, Ai, Ap, gamma) ;
                if ( x [i] >= VSONE )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %25.15e violates upper bound\n",
                                 i, x [i]) ;
                    }
                    two_hit = VSTRUE ;
                    x [i] = VSONE ;
                }
                else dumpj = VSTRUE ;
            }
        }
        else     /* tlinear < 0 */
        {
            /* neg s step */
            s3 = xi/aj ;       /* xi hits lower bound */
            s4 = (1-xj)/ai ;   /* xj hits upper bound */
            if ( s3 < s4 ) /* s is negative and xi hits lower bound */
            {
                pi = -xi ;
                x [i] = VSZERO ;
                Update_Gradient(pi, i, gy, Ai, Ap, gamma) ;
                pj = s3*ai ;
                x [j] += pj ;
                Update_Gradient(pj, j, gy, Ai, Ap, gamma) ;
                if ( x [j] >= VSONE )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %25.15e violates upper bound\n",
                                 j, x [j]) ;
                    }
                    two_hit = VSTRUE ;
                    x [j] = VSONE ;
                }
                else dumpi = VSTRUE ;
            }
            else /* s is negative and xj hits upper bound */
            {
                pj = VSONE - xj ;
                x [j] = VSONE ;
                Update_Gradient(pj, j, gy, Ai, Ap, gamma) ;
                pi = -s4*aj ;
                x [i] += pi ;
                Update_Gradient(pi, i, gy, Ai, Ap, gamma) ;
                if ( x [i] <= VSZERO )
                {
                    if ( PrintLevel )
                    {
                        printf ("x [%i]: %25.15e violates lower bound\n",
                                 i, x [i]) ;
                    }
                    two_hit = VSTRUE ;
                    x [i] = VSZERO ;
                }
                else dumpj = VSTRUE ;
            }
        }

        /* update i and j if two_hit is TRUE */
        if ( two_hit )
        {
            /* find first index > j */
            i = VSP_interior_and_penalty_reduction
                                            (j, x, gx, gy, n, Ai, Ap, gamma) ;
            if ( i == VSEMPTY )
            {
                /* x has no interior components */
                return (EMPTY) ;
            }

            j = VSP_interior_and_penalty_reduction
                                            (i, x, gx, gy, n, Ai, Ap, gamma) ;
            if ( j == VSEMPTY )
            {
                /* done since x has at most one fractional component */
                return (i) ;
            }
            two_hit = VSFALSE ;
        }
        else if ( dumpi )
        {
            i = j ;
            j = VSP_interior_and_penalty_reduction
                                            (j, x, gx, gy, n, Ai, Ap, gamma) ;
            if ( j == VSEMPTY ) return (i) ;
            dumpi = VSFALSE ;
        }
        else /* dumpj */
        {
            j = VSP_interior_and_penalty_reduction
                                            (j, x, gx, gy, n, Ai, Ap, gamma) ;
            if ( j == VSEMPTY ) return (i) ;
            dumpj = VSFALSE ;
        }
    }
}

static VSFLOAT VSP_Obj
(
    VSFLOAT *x,
    VSFLOAT *gx,
    VSFLOAT *y,
    VSFLOAT *a,
    VSINT    n
)
{
    VSINT i ;
    VSFLOAT obj = 0 ;

    for(i = 0 ; i < n ; i++)
    {
        obj += x[i] * gx[i] + y[i] * a[i] ;
    }
    return obj ;
}



/* Compute gz = a - gamma * Aw */
void VSP_Grad
(
    VSFLOAT       *const gz,
    const VSFLOAT *const a,
    const VSFLOAT *const z,
    const VSPINT  *const Ai,
    const VSPLONG *const Ap,
    const VSFLOAT        gamma,
    const VSPINT         n
)
{
    VSPINT col = 0 ;
    VSPLONG p, start, end ;

    /* Initialize Az to all 0's */
    for(col = 0; col < n; col++) gz[col] = a[col] ;

    for (col = 0; col < n; col++)
    {
        VSFLOAT const zcol = -gamma * z[col] ;
        if ( zcol )
        {
            start = Ap[col] ;
            end   = Ap[col + 1] ;
            for (p = start ; p < end ; p++) gz[ Ai[p] ] += zcol ;
            
            /* Take identity into account */
            gz[col] += zcol ;
        }
    }
}

void VSP_Grad_Proj
(
    VSFLOAT  *x,
    VSFLOAT  *gx,
    VSFLOAT  *y,
    VSFLOAT  *gy,
    VSFLOAT  *a,    /* vertex costs */
    VSINT    *Ai,
    VSLINT   *Ap,
    VSINT     u,
    VSINT     n,
    VSFLOAT   gamma,
    int       ngp
)
{
    int it, skip_napheapx, skip_napheapy ;
    VSPINT i ;
    VSFLOAT eps, obj, obj_0, s, t, lambdax, lambday, *X, *Y, *sx, *sy, *gxtemp ;
    int status ;

    NAPdata *napdata ;
    napdata = napheap_setup () ;
    napdata->n = n ;   // problem dimension at current level
    napdata->a = a ;   // linear constraint for napheap problem

#if 0
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
    napdata->blo = 0 ; // lower bound of balance constraint
    napdata->bhi = u ; // upper bound of balance constraint
    eps = n*1.e-8 ;    // when cost does not improve by eps, exit

    X = malloc (n*sizeof (VSPFLOAT)) ; // stores projection in x
    Y = malloc (n*sizeof (VSPFLOAT)) ; // stores projection in y

    // Compute initial objective, obj

    obj = VSP_Obj(x, gx, y, a, n) ;
    lambdax = VSZERO ;
    lambday = VSZERO ;

    /* Perform ngp iterations of the gradient projection algorithm
       The new (x, y) is the solution of
       
              min 0.5*(x-sx)^2 + 0.5*(y-sy)%2
      
                   subject to 0 <= x <= u, 0 <= x <= 1, 0 <= y <= u, 0 <=
    
       where and sx and sy are steps along the gradient starting from the
       iterates x and y: sx = x + s*gx and sy = y + s*gy, s = stepsize */
    sx     = malloc(n * sizeof(*sx)) ;
    sy     = malloc(n * sizeof(*sy)) ;
    gxtemp = malloc(n * sizeof(*gxtemp)) ;

    for(it = 1 ; it <= ngp ; it++)
    {
        obj_0 = obj ; // save prior objective value
        s = vs_max(gx, n) ;
        t = vs_max(gy, n) ;
        s = MAX(s, t) ;

        // Scaling factor for gradients
        s = 1 / s ;

        /* Project gradient step onto feasible set. The projection decomposes
           into two separate projections: sx is projected on the x-feasible
           set and sy is projected on the y-feasible set.  */
        obj = -NAPINF ; /* initialize obj to a lower bound on objective*/
        int subit = 0 ;
        while(obj + eps < obj_0 )
        {
            subit++ ; 
            // Scale gradients by 0.5
            if ( subit == 1 )
            {
                vs_scale(gx, VSHALF*s, n) ;
                vs_scale(gy, VSHALF*s, n) ;
            }
            else
            {
                vs_scale(gx, VSHALF, n) ;
                vs_scale(gy, VSHALF, n) ;
            }

            // Compute new point (sx, sy)
            vs_vecadd(sx, x, gx, n) ;
            vs_vecadd(sy, y, gy, n) ;

            // Project sx to feasible set, see if lambda = 0 works
            skip_napheapx = VSFALSE ;
            if ( !lambdax )
            {
                t = VSZERO ;
                for (i = 0; i < n; i++)
                {
                    if ( sx [i] > VSONE )
                    {
                        sx [i] = VSONE ;
                        t += a [i] ;
                    }
                    else if ( sx [i] < VSZERO )
                    {
                        sx [i] = VSZERO ;
                    }
                    else
                    {
                        t += a [i]*sx [i] ;
                    }
                }
                if ( t <= u )
                {
                    skip_napheapx = VSTRUE ;
                }
                else
                {
                    vs_vecadd(sx, x, gx, n) ; /* sx restored to correct value*/
                }
            }
            if ( !skip_napheapx )
            {
                napdata->c = sx ;
                napdata->x = X ;

                napdata->lambda = lambdax ;
                status = napheap(napdata) ;
                lambdax = napdata->lambda ;
#if 0
                //printf ("lambdax: %e\n", lambdax) ;
                //if( status ) napheap_print_status (napdata) ;
                VSFLOAT err, r ;            
                err = napheap_check (&r, &s, &t, napdata) ;
                printf ("errx: %e errb: %e errB: %e errc: %e\n", err, r, s, t) ;
                if ( err > 1.e-6 )
                {
                    printf ("NAPHEAP error: %e\n", err) ;
                    exit(1) ;
                }
#endif
            }

            // Project sy to feasible set, see if lambda = 0 works
            skip_napheapy = VSFALSE ;
            if ( !lambday )
            {
                t = VSZERO ;
                for (i = 0; i < n; i++)
                {
                    if ( sy [i] > VSONE )
                    {
                        sy [i] = VSONE ;
                        t += a [i] ;
                    }
                    else if ( sy [i] < VSZERO )
                    {
                        sy [i] = VSZERO ;
                    }
                    else
                    {
                        t += a [i]*sy [i] ;
                    }
                }
                if ( t <= u )
                {
                    skip_napheapy = VSTRUE ;
                }
                else
                {
                    vs_vecadd(sy, y, gy, n) ; /* sy restored to correct value*/
                }
            }
            if ( !skip_napheapy )
            {
                napdata->c = sy ;
                napdata->x = Y ;
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
                printf ("erry: %e errb: %e errB: %e errc: %e\n", err, r, s, t) ;
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
                VSP_Grad(gxtemp, a, sy, Ai, Ap, gamma, n) ;
                if ( skip_napheapx )
                {
                    // solution stored in sx
                    obj = VSP_Obj(sx, gxtemp, sy, a, n) ;
                }
                else
                {
                    // solution stored in X
                    obj = VSP_Obj(X, gxtemp, sy, a, n) ;
                }
            }
            else
            {
                // solution stored in Y
                VSP_Grad(gxtemp, a, Y, Ai, Ap, gamma, n) ;
                if ( skip_napheapx )
                {
                    // solution stored in sx
                    obj = VSP_Obj(sx, gxtemp, Y, a, n) ;
                }
                else
                {
                    // solution stored in X
                    obj = VSP_Obj(X, gxtemp, Y, a, n) ;
                }
            }

            // printf("it: %i obj: %f    obj0: %f\n", it, obj, obj_0) ;
        }
        

        /* Update: x <- X and y <- Y */
        if ( skip_napheapx )
        {
            // solution stored in sx
            vs_copyx(x, sx, n) ;
        }
        else
        {
            // solution stored in X
            vs_copyx(x, X, n) ;
        }
        if ( skip_napheapy )
        {
            // solution stored in sy
            vs_copyx(y, sy, n) ;
        }
        else
        {
            // solution stored in Y
            vs_copyx(y, Y, n) ;
        }

        // grad at x stored above in gxtemp
        vs_copyx(gx, gxtemp, n) ;
        
        // evaluate gradient at y
        VSP_Grad(gy, a, x, Ai, Ap, gamma, n) ;

        if ( obj < obj_0 + eps ) break ; // |obj - obj_0| < eps => stationary

    }

    napheap_terminate(&napdata) ;

    /* Cleanup */
    free(sx) ;
    free(sy) ;
    free(X) ;
    free(Y) ;
    free(gxtemp) ;
}


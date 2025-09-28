#include "shared_api.h"
#include "ds_api.h"
#include <string.h>
#include "Stopwatches.h"
#include "shared_types.h"


/******************* VSP_stats ************************************************/
VSP_stats *VSP_construct_stats()
{
    /*--- Generic Info -------------------------------------------------------*/
    VSP_stats *stats = malloc( sizeof(*stats) ) ;

    /*--- Resolve Level ------------------------------------------------------*/
    stats->resolve_tried  = VSP_construct_dyn_stk() ;
    stats->resolve_better = VSP_construct_dyn_stk() ;
    
    stats->objIn  = VSP_construct_dyn_stk() ;
    stats->objOut = VSP_construct_dyn_stk() ;
    stats->sepIn  = VSP_construct_dyn_stk() ;
    stats->sepOut = VSP_construct_dyn_stk() ;
    stats->numPen = VSP_construct_dyn_stk() ;

    return stats ;
}

void VSP_destruct_stats(VSP_stats *stats)
{
    VSP_destruct_dyn_stk(stats->resolve_tried) ;
    VSP_destruct_dyn_stk(stats->resolve_better) ;
    VSP_destruct_dyn_stk(stats->objIn) ;
    VSP_destruct_dyn_stk(stats->objOut) ;
    VSP_destruct_dyn_stk(stats->sepIn) ;
    VSP_destruct_dyn_stk(stats->sepOut) ;
    VSP_destruct_dyn_stk(stats->numPen) ;
    free(stats) ;
}


/**************** VSP_solve_result ********************************************/
void VSP_destruct_solve_result(VSP_solve_result *result)
{
    free(result->x) ;
    free(result->y) ;
    free(result) ;

    return ;
}

/**************** VSP_point ***************************************************/
VSP_point* VSP_construct_point(VSPINT n) {
    VSP_point* point = malloc( sizeof(VSP_point) ) ;
    point->int_part  = calloc( n, sizeof(VSPINT) ) ; // calloc to initialize 0

    /* Initialize default values */
    point->ind_frac = EMPTY ;
    point->val_frac = 0 ;
    point->n = n ;

    return point ;
}

void VSP_destruct_point(VSP_point* point)
{
    #ifdef VSP_DEBUG
        static long long unsigned count = 0 ;
        count ++ ;
    #endif

    free(point->int_part) ;
    free(point) ;
    return ;
}



/**************** VSP_Az ******************************************************/
VSP_Az* VSP_construct_Az(VSPINT n) {
    VSP_Az* Az     = malloc( sizeof(VSP_Az) ) ;
    Az->int_part   = calloc(n, sizeof(VSPINT)) ; // calloc to 0
    Az->col_frac_z = calloc(n, sizeof(VSPINT)) ;    // calloc to 0

    /* Initialize default values */
    Az->ind_frac_z = EMPTY ;
    Az->val_frac_z = 0 ;
    Az->one = 1 ;
    Az->n = n ;

    return Az ;
}

void VSP_destruct_Az(VSP_Az* Az) {
    free(Az->int_part) ;
    free(Az->col_frac_z) ;
    free(Az) ;
}

/************** VSP_flips_dec *************************************************/
VSP_flips_dec *VSP_construct_flips_dec(VSP_flips_dec_type decType)
{
    VSP_flips_dec *ret = calloc(1, sizeof( *ret) ) ;
    ret->decType = decType ;
    ret->levels = NULL ;
    ret->decBoundary = EMPTY ;

    return ret ;
}

void VSP_destruct_flips_dec(VSP_flips_dec *dec)
{
    if(dec->levels != NULL)
    { 
        free(dec->levels) ;
    }
    free(dec) ;

    return ;
}



/***************** VSP_state **************************************************/

VSP_state* VSP_construct_state(VSPINT n)
{
    VSP_state* Z = malloc(sizeof(VSP_state)) ;
    Z->point     = VSP_construct_point(n) ;
    Z->Az        = VSP_construct_Az(n) ;

    /* Set default values */
    Z->weight    = 0 ;

    return Z ;
}

void VSP_destruct_state(VSP_state* Z)
{
    VSP_destruct_point(Z->point) ;
    VSP_destruct_Az(Z->Az) ;
    free(Z) ;
    return ;
}


/***************** VSP_coarsen_dec ********************************************/
/* Does NOT construct levels struct */
VSP_coarsen_dec *VSP_construct_coarsen_dec(VSP_coarsen_dec_type decType)
{
    VSP_coarsen_dec *ret = calloc( 1, sizeof(*ret) ) ;
    ret->decType = decType ;
    ret->levels  = NULL ;
    
    return ret ;
}

void VSP_destruct_coarsen_dec(VSP_coarsen_dec *dec)
{
    if(dec->levels != NULL)
    {
        VSP_destruct_levels(dec->levels) ;
    }
    free(dec) ;
}


/***************** VSP_options ************************************************/
VSP_options* VSP_construct_options()
{
    VSP_options* options = malloc(sizeof(VSP_options));
    options->restart_decision     = NULL ;
    return options;
}


void VSP_destruct_options(VSP_options* options)
{
    if(options->restart_decision != NULL)
    {
        VSP_destruct_restart_decision(options->restart_decision) ;
    }

    if(options->coarsenDec != NULL)
    {
        VSP_destruct_coarsen_dec(options->coarsenDec) ;
    }

    /* Free the struct itself */
    free(options);
}


// Allocates memory but does NOT set it to anything
// The max vertex cost, min vertex cost, and total vertex cost are set to 0
VSP_Graph *Construct_Graph(VSPINT n, VSPLONG nz)
{
    VSP_Graph *Graph = calloc(1, sizeof(*Graph)) ;

    Graph->n         = n  ;
    Graph->nz        = nz ;

    Graph->Ae        = malloc(  nz     * sizeof(VSPLONG) ) ;
    Graph->Ai        = malloc(  nz     * sizeof(VSPLONG) ) ;
    Graph->Ap        = malloc( (n + 1) * sizeof(VSPINT)  ) ;

    Graph->c         = malloc(  n      * sizeof(VSPINT)  ) ;
    return Graph ;
}

void Destroy_Graph(VSP_Graph *Graph)
{
    free(Graph->Ae) ;
    free(Graph->Ai) ;
    free(Graph->Ap) ;
    free(Graph->c ) ;
    free(Graph) ;
}


/* *************** Begin VSP_problem **************************************** */
VSP_problem* VSP_construct_problem
(
    VSPINT       num_vertices,
    VSPINT                 nz,
    VSP_options      *options
)
{
    /* Create the problem */
    VSP_problem* problem = calloc(1, sizeof(VSP_problem));

    /* Construct the name */
    problem->name = malloc( 350 * sizeof(char) ) ;
    if(problem->name == NULL)
    {
        printf("\n\nOUT OF MEMORY FOR NAME\n\n") ;
        exit(1) ;
    }
    problem->name[0] = '\0' ;

    /* Construct adjacency matrix */
    problem->Graph     = Construct_Graph(num_vertices, nz) ;

    /* Construct matching arrays. Need to use calloc */
    problem->matching = calloc( num_vertices, sizeof(VSPINT));
    problem->matchmap = calloc( num_vertices, sizeof(VSPINT));
    problem->invmatchmap = calloc( num_vertices, sizeof(VSPINT));
    problem->matchtype = calloc( num_vertices, sizeof(VSP_match_type));

    /* Initialize scalars */
    problem->Graph->n  = num_vertices;
    problem->Graph->nz = nz;
    problem->cn = 0;
    problem->cnz = 0;
    problem->singleton = EMPTY ;
    problem->parent_problem = NULL ;
    problem->child_problem  = NULL ;
    problem->level = 0 ;
    problem->coarsest_level_num = 0 ;
    problem->ratioMaxCostNodes = 0 ;

    problem->maxDegree = 0 ;

    problem->resolveTried   = false ;
    problem->resolveSuccess = false ;

    problem->stats = VSP_construct_stats() ;

    // Construct Stopwatches
    problem->Stopwatches = Stopwatches_Construct(30, problem) ;

    /* The graph is assumed to be uncompressed */
    problem->Graph         = problem->Graph ;

    return problem;
}

void VSP_destruct_problem(VSP_problem* problem)
{
    /* Free the graph data */
    free(problem->name) ;

    Destroy_Graph(problem->Graph) ;

    /* Free matching arrays */
    free(problem->matching);
    free(problem->matchmap);
    free(problem->invmatchmap);
    free(problem->matchtype);


    VSP_destruct_stats(problem->stats) ;
    Stopwatches_Destroy(problem->Stopwatches) ;

    /* Free the struct itself */
    /* Note: we are not freeing cproblem, since each problem frees itself*/
    /* via the following line */
    free(problem);
}


/* ************* Explicit Matrix ******************************************* **/
VSPINT** VSP_construct_explicit_matrix
(
    VSPINT *values,
    VSPINT     *Ai,
    UNUSED VSPINT      nz,
    VSPINT     *Ap,
    VSPINT       n
)
{
    static long long unsigned count = 0 ;
    count += 1 ;

    VSPINT i, j, p, end ;

    /* Make pointers for the rows of A */
    VSPINT** A = malloc( n * sizeof(VSPINT*) ) ;
    if ( A == NULL ) {
        VSP_error("Memory allocation failed", DBG_INFO) ;
        exit(1) ;
    }

    /* Allocate space for each row */
    for (i = 0 ; i < n ; i++) {
        A [i] = malloc( n * sizeof(VSPINT) ) ;
        if ( A [i] == NULL ) {
            VSP_error("Memory allocation failed", DBG_INFO) ;
            exit(1) ;
        }
    }

    /* Zero out the matrix */
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) A [i][j] = 0 ;
    }

    /* Fill out the matrix with 1's from Ai and Ap */
    for (j = 0 ; j < n ; j++) {
        end = Ap [ j + 1 ] ;
        for (p = Ap[j] ; p < end ; p++) {
            i = Ai [p] ;
            if (values == NULL) A [i][j] = 1;
            else A [i][j] = values [p] ;
        }
    }

    /* Check that the matrix is symmetric */
    int is_symmetric = TRUE ;
    for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
            if ( A[i][j] != A [j][i] ) {
                is_symmetric = FALSE ;
                printf("\nA[%i][%i] = %i   while   A[%i][%i] = %i", i, j, A[i][j], j, i, A[j][i]) ;
                fflush(stdout) ;
            }
        }
    }

    if (is_symmetric == FALSE)
    {
        if (values == NULL)
        {
            VSP_error("Adjacency matrix not symmetric", DBG_INFO) ;
        }
        else
        {
            VSP_error("Edge weights matrix not symmetric", DBG_INFO) ;
        }
    }

    return A ;
}

void VSP_destruct_explicit_matrix( VSPINT** A, VSPINT n ) {
    VSPINT i ;

    for ( i = 0 ; i < n ; i++ ) free (A [i]) ;
    free(A) ;

    return ;
}
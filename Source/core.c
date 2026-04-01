#include "core.h"



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


/* *************** Graph **************************************************** */
VSP_Graph *Construct_Graph(VSPINT n, VSPLONG nz, VSPLONG *Ap, VSPINT *Ai, VSPINT*Ae)
{
    VSP_Graph *Graph = calloc(1, sizeof(*Graph)) ;

    Graph->n         = n  ;
    Graph->nz        = nz ;

    if(Ae == NULL)
    {
        Graph->Ae = malloc(  nz     * sizeof(VSPINT) ) ;
    }
    else Graph->Ae = Ae ;

    if(Ai == NULL)
    {
        Graph->Ai = malloc(  nz     * sizeof(VSPINT) ) ;
    }
    else Graph->Ai = Ai ;

    if(Ap == NULL)
    {
        Graph->Ap = malloc( (n + 1) * sizeof(VSPLONG)  ) ;
    }
    else Graph->Ap = Ap ;

    Graph->c = malloc(  n      * sizeof(VSPINT)  ) ;
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


/* *************** VSP_problem ********************************************** */
VSP_problem* VSP_construct_problem
(
    VSPINT       num_vertices,
    VSPLONG                nz,
    VSPINT                *Ai, // If NULL, then this is allocated
    VSPLONG               *Ap, // If NULL, then this is allocated
    VSPINT                *Ae, // If NULL, then this is allocated
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
    problem->Graph     = Construct_Graph(num_vertices, nz, Ap, Ai, Ae) ;

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
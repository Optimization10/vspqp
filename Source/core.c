#include "core.h"
#include <stdio.h>


/**************** Point ***************************************************/
Point* Construct_Point(VSPINT n) {
    Point* point = malloc( sizeof(Point) ) ;
    point->int_part  = calloc( n, sizeof(VSPINT) ) ; // calloc to initialize 0

    /* Initialize default values */
    point->ind_frac = EMPTY ;
    point->val_frac = 0 ;
    point->n = n ;

    return point ;
}

void Destroy_Point(Point* point)
{
    #ifdef DEBUG
        static long long unsigned count = 0 ;
        count ++ ;
    #endif

    free(point->int_part) ;
    free(point) ;
    return ;
}



/**************** Az ******************************************************/
Az* Construct_Az(VSPINT n) {
    Az* Az     = malloc( sizeof(*Az) ) ;
    Az->int_part   = calloc(n, sizeof(VSPINT)) ; // calloc to 0
    Az->col_frac_z = calloc(n, sizeof(VSPINT)) ;    // calloc to 0

    /* Initialize default values */
    Az->ind_frac_z = EMPTY ;
    Az->val_frac_z = 0 ;
    Az->one = 1 ;
    Az->n = n ;

    return Az ;
}

void Destroy_Az(Az* Az)
{
    free(Az->int_part) ;
    free(Az->col_frac_z) ;
    free(Az) ;
}

/***************** State **************************************************/

State* Construct_State(VSPINT n)
{
    State* Z = malloc(sizeof(State)) ;
    Z->point     = Construct_Point(n) ;
    Z->Az        = Construct_Az(n) ;

    /* Set default values */
    Z->weight    = 0 ;

    return Z ;
}

void Destroy_State(State* Z)
{
    Destroy_Point(Z->point) ;
    Destroy_Az(Z->Az) ;
    free(Z) ;
    return ;
}


/* *************** Graph **************************************************** */
Graph *Construct_Graph(VSPINT n, VSPLONG nz, VSPLONG *Ap, VSPINT *Ai, VSPINT*Ae)
{
    Graph *Graph = calloc(1, sizeof(*Graph)) ;

    Graph->n         = n  ;
    Graph->nz        = nz ;

    if(Ae == NULL)
    {
        Graph->Ae = malloc(  nz     * sizeof(VSPINT) ) ;
        Graph->owns_Ae = true ;
    }
    else
    {
        Graph->Ae = Ae ;
        Graph->owns_Ae = false ;
    }

    if(Ai == NULL)
    {
        Graph->Ai = malloc(  nz     * sizeof(VSPINT) ) ;
        Graph->owns_Ai = true ;
    }
    else
    {
        Graph->Ai = Ai ;
        Graph->owns_Ai = false ;
    }

    if(Ap == NULL)
    {
        Graph->owns_Ap = true ;
        Graph->Ap = malloc( (n + 1) * sizeof(VSPLONG)  ) ;
    }
    else
    {
        Graph->owns_Ap = false ;
        Graph->Ap = Ap ;
    }

    Graph->c = malloc(  n      * sizeof(VSPINT)  ) ;
    return Graph ;
}

void Destroy_Graph(Graph *Graph)
{
    if(Graph->owns_Ae) free(Graph->Ae) ;
    if(Graph->owns_Ai) free(Graph->Ai) ;
    if(Graph->owns_Ap) free(Graph->Ap) ;
    free(Graph->c ) ;
    free(Graph) ;
}

/* *************** Matching ************************************************* */
/* Note: matching is not really a core module. However, it is here because
   because the problem destroyer should call the matching destroyer, and if
   I make matching a part of the coarsening module, I would have to include the
   coarsening api into core .c, which I don't want to do since */
Matching *Matching_Construct(VSPINT n)
{
    /* Construct matching arrays. Need to use calloc */
    Matching *matching    = malloc(sizeof(*matching)) ;
    matching->match       = calloc( n, sizeof(VSPINT)) ;
    matching->matchmap    = calloc( n, sizeof(VSPINT)) ;
    matching->invmatchmap = calloc( n, sizeof(VSPINT)) ;
    matching->matchtype   = calloc( n, sizeof(match_type)) ;
    matching->cn          = 0 ;
    return matching ;
}

void Matching_Destroy(Matching *matching)
{
    free(matching->match) ;
    free(matching->matchmap) ;
    free(matching->invmatchmap) ;
    free(matching->matchtype) ;
    free(matching) ;
}


/* *************** Problem ************************************************** */
Problem* Construct_Problem
(
    VSPINT       num_vertices,
    VSPLONG                nz,
    VSPINT                *Ai, // If NULL, then this is allocated
    VSPLONG               *Ap, // If NULL, then this is allocated
    VSPINT                *Ae, // If NULL, then this is allocated
    Options          *options
)
{
    /* Create the problem */
    Problem* problem = calloc(1, sizeof(Problem));

    /* Construct adjacency matrix */
    problem->graph = Construct_Graph(num_vertices, nz, Ap, Ai, Ae) ;

    /* Initialize the matching to NULL */
    problem->matching = Matching_Construct(num_vertices) ;

    /* Initialize scalars */
    problem->graph->n  = num_vertices;
    problem->graph->nz = nz;
    problem->parent_problem = NULL ;
    problem->child_problem  = NULL ;
    problem->level = 0 ;
    problem->coarsest_level_num = 0 ;

    return problem;
}

// Npte: does NOT free matching
void Destroy_Problem(Problem* problem)
{
    Destroy_Graph(problem->graph) ;
    Matching_Destroy(problem->matching) ;
    free(problem);
}


void Error
( 
    char* msg,
    char* file_name,
    VSPINT line_num,
    unsigned long long count
)
{
    #ifndef DEBUG
        printf("\n\nWARNING! DOING DEBUGGING CHECKS WHILE NOT IN DEBUG MODE") ;
    #endif

    #ifdef DEBUG
        static long long unsigned error_count = 0 ;
        error_count++ ;


        FILE* file = fopen ( "output/Error_log.txt", "a" ) ;
        if(error_count <= 10000)
        {
            fprintf(file, "\n\n ERROR: %s", msg) ;
            if(G_Current_Problem_Name != NULL)
            {
                fprintf(file, "\n Prob name   : %s", G_Current_Problem_Name) ;
            }
            else
            {
                fprintf(file, "\n Prob name   : Not Provided") ;
            }
            fprintf(file, "\n File name   : %s", file_name) ;
            fprintf(file, "\n Line number : %lld", (long long) line_num) ;
            fprintf(file, "\n Count : %llu", count) ;
            fclose(file) ;
        }
        else
        {
            fprintf(file, "\n\n: Error log too large, terminating program.\n") ;
            fclose(file) ;
            printf("\n\nERROR: Error log too large, terminating program.\n\n") ;
            exit(1) ;
        }
    #endif
    return ;
}
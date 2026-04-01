#ifndef CORE_H
#define CORE_H

#include <time.h>
#include <stdbool.h>
// #include "VSPQP.h"


/* ========================================================================== */
// Core Types
/* ========================================================================== */

typedef enum VSP_match_type_enum
{
    matchtype_unmatched = 0,
    matchtype_orphan,
    matchtype_standard,
    matchtype_brotherly,
    matchtype_community,

    // Matchtypes for even coarsening
    matchtype_selfMatch,
    matchtype_singleton,
    matchtype_HEM,
    matchtype_3way
} VSP_match_type ;


typedef struct VSP_problem_struct VSP_problem;

typedef struct
{
    VSPINT  n ;
    VSPLONG nz ;

    VSPINT  *Ai ;
    VSPLONG *Ap ;
    VSPINT  *Ae ;

    VSPINT *c ;
    VSPINT min_vtx_cost ;
    VSPINT max_vtx_cost ;
    VSPINT total_vtx_cost ;
} VSP_Graph ;

typedef struct Stopwatch
{
    char           *name ;
    struct timespec start ;
    struct timespec end ;
    bool            inProgress ;
    long double     accumTime ;
    long long   counter ; // Number of times the sw has been start/stopped
} Stopwatch ;

typedef struct Stopwatches
{
    size_t        numWatches ;
    size_t        maxNumWatches ;
    Stopwatch   **watches ;
    VSP_problem  *problem ; // Associated with this collection of stopwatches
} Stopwatches ;

typedef struct VSP_problem_struct
{
    // This is the graph which is being solved
    VSP_Graph             *Graph ;

    char                   *name ;  /* name of the graph                      */

    /* Further problem info */
    VSPINT                     u ;  /* upper balance constraint               */

    /* Info for matching */
    VSPINT             *matching ;
    VSPINT             *matchmap ;
    VSPINT          *invmatchmap ;
    VSP_match_type    *matchtype ;
    VSPINT             singleton ;

    /* Coarser graph size */
    VSPINT                    cn ;
    VSPLONG                   cnz ;

    /* Link to parent problem. NULL if this is finest graph */
    struct VSP_problem_struct* parent_problem ;
    bool was_coarsened_w_comm ;

    /* Link to child problem. NULL if this is the coarsest graph */
    struct VSP_problem_struct* child_problem ;

    /* The number of coarsenings which were done to reach this Graph-> The
       finest problem is at level 0 */
    VSPINT level ;
    VSPINT coarsest_level_num ; /* The level of the coarsest problem */

    /* Re-solve */
    bool resolveTried ;
    bool resolveSuccess ;

    VSPINT ubFinest ;

    /* The ratio of nodes which have maximal cost*/
    VSPFLOAT ratioMaxCostNodes ;

    /* The maximum degree in the graph */
    VSPINT maxDegree ;

    VSP_stats *stats ;

    Stopwatches *Stopwatches ;

} VSP_problem ;

/* ========================================================================== */
// Constructors and Destructors for Core Types
/* ========================================================================== */

VSP_problem* VSP_construct_problem
(
    VSPINT       num_vertices,
    VSPLONG                nz,
    VSPINT                *Ai, // If NULL, then this is allocated
    VSPLONG               *Ap, // If NULL, then this is allocated
    VSPINT                *Ae, // If NULL, then this is allocated
    VSP_options      *options
) ;

void VSP_destruct_problem(VSP_problem* problem) ;



VSP_Graph *Construct_Graph(VSPINT n, VSPLONG nz, VSPLONG *Ap, VSPINT *Ai, VSPINT*Ae) ;
void Destroy_Graph(VSP_Graph *Graph) ;

VSP_Az* VSP_construct_Az(VSPINT n) ;
void VSP_destruct_Az(VSP_Az* Az) ;

VSP_point* VSP_construct_point(VSPINT n ) ;
void VSP_destruct_point(VSP_point* point) ;

VSP_state* VSP_construct_state(VSPINT n);
void VSP_destruct_state(VSP_state* state);

#endif
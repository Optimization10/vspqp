#ifndef CORE_H
#define CORE_H

#include <stdbool.h>
#include <math.h>
#include <stdlib.h>  
#include <stdio.h>
#include "VSPQP.h"


/* ========================================================================== */
/*  Core Macros                                                               */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/*  Constants                                                                 */
/* -------------------------------------------------------------------------- */

#define EMPTY  -1
#define PLUS    1
#define MINUS  -1
#define UP      1
#define DOWN   -1

#define TRUE    1
#define FALSE   0

#define VSP_ZERO    1E-6
#define ZERO_PLUS   1E-6
#define ZERO_MINUS -(ZERO_PLUS)

#define VAR_X 0
#define VAR_Y 1

/* -------------------------------------------------------------------------- */
/*  Utility                                                                   */
/* -------------------------------------------------------------------------- */

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#define IN_RANGE(num, low, high) ( ((low) <= (num)) && ((num) < (high)) )

#define IFF(A, B) ( ((A) + (B) == 0) || ((A) + (B) == 2) )

#define CAST(arg, type)  ((type)(arg))
#define PRI(A)           ((long long int)(A))

#define UNUSED  __attribute__((unused))

/* Returns true iff a/b < c/d  (and <= variant) */
#define CMP_FRAC(a, b, c, d)    ((a) * (d) <  (c) * (b))
#define CMPEQ_FRAC(a, b, c, d)  ((a) * (d) <= (c) * (b))

#define IMBALANCE(wgtz, wgtw) \
({                             \
    typeof(wgtz) _diff = (wgtz) - (wgtw); \
    _diff < 0 ? -_diff : _diff;            \
})

#define SWAP(a, b, type) do {                  \
    type _swap_tmp_do_not_use_ = (a);          \
    (a) = (b);                                 \
    (b) = _swap_tmp_do_not_use_;               \
} while (0)

#define ABS_VAL(a) ( ((a) >= 0 ) ? ( (a) ) : ( -(a) ) )

/* -------------------------------------------------------------------------- */
/*  Debug                                                                     */
/* -------------------------------------------------------------------------- */

#define ERR_MSG_LEN 300

#define DBG_INFO  __FILE__, __LINE__, count
#define ERR(msg)  Error(msg, __FILE__, __LINE__, count)

#ifdef DEBUG
    #define IFDBG(code)               (code)
    #define SET_FUNCTION_COUNTER      static long long unsigned count = 0; count++
    #define VSP_ASSERT(cond) do {                              \
        if (!(cond)) Error("Failed assert.", DBG_INFO);        \
    } while (0)
#else
    #define IFDBG(code)
    #define SET_FUNCTION_COUNTER
    #define VSP_ASSERT(cond)
#endif

#if defined(DEBUG) || defined(DEBUG_COARSENING)
    #define WARN_IF_NOT_DEBUG_MODE
#else
    #define WARN_IF_NOT_DEBUG_MODE \
        printf("\n\nWARNING: Calling debug function outside debug mode. " \
               "File: %s  Line: %lld\n\n", __FILE__, (long long)__LINE__)
#endif


/* ========================================================================== */
/*  Core Types                                                                */
/* ========================================================================== */

typedef struct
{
    VSPINT  n ;
    VSPLONG nz ;

    VSPINT  *Ai ;
    VSPLONG *Ap ;
    VSPINT  *Ae ;

    VSPINT   *c ;
    VSPINT    max_deg ;
    VSPINT    min_vtx_cost ;
    VSPINT    max_vtx_cost ;
    VSPINT    total_vtx_cost ;
    VSPFLOAT  ratio_max_cost ; // ratio of vertices with max cost

    VSPINT singletons_cost ; // EMPTY if this was not computed

    // If true, then memory was allocated by the graph constructor
    // which needs to be freed by the destroyer
    bool owns_Ai ;
    bool owns_Ap ;
} Graph ;

typedef enum match_type_enum
{
    matchtype_unmatched = 0,
    matchtype_orphan,
    matchtype_standard,
    matchtype_brotherly,
    matchtype_community,

    /* Matchtypes for even coarsening */
    matchtype_selfMatch,
    matchtype_singleton,
    matchtype_HEM,
    matchtype_3way
} match_type ;

typedef struct
{
    VSPINT           *match ;
    VSPINT        *matchmap ;
    VSPINT     *invmatchmap ;
    match_type   *matchtype ;
    VSPINT               cn ; // Number of vertices in coarse graph
} Matching ;

typedef struct stats
{
    // EMPTY if no restart tried, 0 if tried and failed,1 if tried and succeeded
    int      restart_success ;
    VSPFLOAT restart_obj_ratio ; // new obj / old obj
} Stats ;

typedef struct Problem_struct
{
    Graph *graph ;  /* the graph being solved */
    VSPINT     u ;  /* upper balance constraint */

    /* Matching info */
    Matching *matching ;

    /* Pointers to parent and child problem */
    struct Problem_struct *parent_problem ;
    struct Problem_struct *child_problem ;

    /* Number of coarsenings done to reach this graph Finest is level 0*/
    VSPINT level ;
    VSPINT coarsest_level_num ; /* the level of the coarsest problem */

    /* Vertices in the input graph, removed hubs included. Hubs are expected to
       stay in the separator, so u at every level is a ratio of this count. */
    VSPINT n_original ;

    Stats stats ;

} Problem ;

/* Store an almost binary point */
typedef struct
{
    VSPINT    *int_part ; /* integer part                                     */
    VSPINT     ind_frac ; /* index of fractional entry. EMPTY if none.        */
    VSPFLOAT   val_frac ; /* value of fractional entry. No guarantees
                             if ind_frac is EMPTY.                            */
    VSPINT            n ; /* number of entries in the point                   */
} Point ;

/* Store Az where A is a binary matrix and z has at most one fractional entry.
   Entry i in Az is
   int_part[i] + val_frac_z * ( col_frac_z[i] == one ) * (ind_frac_z != EMPTY)*/
typedef struct
{
    VSPINT     *int_part ; /* A times the integer part of z                   */
    VSPINT    ind_frac_z ; /* index of fractional entry in z. EMPTY if none.  */
    VSPFLOAT  val_frac_z ; /* value of fractional entry in z. No guarantees
                              if ind_frac_z is EMPTY.                         */
    VSPINT   *col_frac_z ; /* column of A corresponding to frac index in z.
                              No guarantees if ind_frac_z is EMPTY.           */
    VSPINT           one ; /* largest integer in col_frac_z. Represents the
                              1's in that column. All smaller integers
                              represent zeros.                                */
    VSPINT             n ; /* number of entries in Az                         */
} Az ;

/* Container for all the information pertaining to x or y */
typedef struct
{
    Point             *point ; /* x or y                                      */
    Az                   *Az ; /* Ax or Ay                                    */
    VSPINT            weight ; /* c'x or c'y                                  */
} State ;
/* ========================================================================== */
// Constructors and Destructors for Core Types
/* ========================================================================== */

Problem* Construct_Problem
(
    VSPINT       num_vertices,
    VSPLONG                nz,
    VSPINT                *Ai, // If NULL, then this is allocated
    VSPLONG               *Ap, // If NULL, then this is allocated
    bool             alloc_Ae,
    Options          *options
) ;

void Destroy_Problem(Problem* problem) ;



Graph *Construct_Graph
(
    VSPINT n, VSPLONG nz, VSPLONG *Ap, VSPINT *Ai, bool alloc_Ae
) ;

void Destroy_Graph(Graph *Graph) ;

Az* Construct_Az(VSPINT n) ;
void Destroy_Az(Az* Az) ;

Point* Construct_Point(VSPINT n ) ;
void Destroy_Point(Point* point) ;

State* Construct_State(VSPINT n);
void Destroy_State(State* state);


void Error
( 
    char*              msg,
    char*              file_name,
    VSPINT             line_num,
    unsigned long long count
) ;

/* ========================================================================== */
// Global Variables
/* ========================================================================== */
extern const char *G_Current_Problem_Name ;


#endif
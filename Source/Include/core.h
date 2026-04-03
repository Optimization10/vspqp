#ifndef CORE_H
#define CORE_H

#include <time.h>
#include <stdbool.h>
#include "VSPQP.h"


/* ========================================================================== */
// Core Macros
/* ========================================================================== */
#define EMPTY -1

    #define TIMER_INIT \
    struct timespec hidden_timespec_struct_start ; \
    struct timespec hidden_timespec_struct_end ;

    #define TIMER_START \
    struct timespec hidden_timespec_struct_start ; \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_start)

    #define TIMER_ELAPSED \
    (VSP_compute_elapsed_time(hidden_timespec_struct_start, hidden_timespec_struct_end))
    
    #define TIMER_END \
    struct timespec hidden_timespec_struct_end ;   \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_end)

    #define PERFORMANCE_UPDATE_TIME(PREFIX) \
    // do { \
    //     G_Performance_Struct.PREFIX##_Count++; \
    //     G_Performance_Struct.PREFIX##_Seconds += \
    //     (   (hidden_timespec_struct_end.tv_sec - hidden_timespec_struct_start.tv_sec) \
    //     + (hidden_timespec_struct_end.tv_nsec - hidden_timespec_struct_start.tv_nsec) / 1e9L ) ;\
    // } while(0)

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


typedef struct Problem_struct Problem;

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





typedef struct Problem_struct
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
    struct Problem_struct* parent_problem ;
    bool was_coarsened_w_comm ;

    /* Link to child problem. NULL if this is the coarsest graph */
    struct Problem_struct* child_problem ;

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

} Problem ;

/* Store an almost binary point.                                              */
typedef struct {
    VSPINT    *int_part ; /* Integer part*/
    VSPINT     ind_frac ; /* Index of fractional entry. EMPTY if none.        */
    VSPFLOAT   val_frac ; /* Value of fractional entry. No guarantees
                             if ind_frac is EMPTY.                            */
    VSPINT            n ; /* Number of entries in the point                   */
} Point ;


/* Store Az where A is a binary matrix and z has at most one fractional entry
   Entry i in Az is
   int_part[i] + val_frac_z * ( col_frac_z[i] == one ) * (ind_frac_z != EMPTY)*/
typedef struct {
    VSPINT     *int_part ; /* A times the integer part of z                   */
    VSPINT    ind_frac_z ; /* Index of fractional entry in z. EMPTY if none.  */
    VSPFLOAT  val_frac_z ; /* Value of fractional entry in z. No guarantees if
                              if ind_frac_z is EMPTY.                       */
    VSPINT   *col_frac_z ; /* Column of A corresponding to frac index in z. No
                              guarantees if ind_frac_z is EMPTY.              */
    VSPINT           one ;    /* Largest integer in col_frac_v. Represents the
                              1'S in that column. All smaller integers represent
                              zeros.                                          */
    VSPINT             n ; /* Number of entries in Az.                        */
} VSP_Az ;

/* Container for all the information pertaining to x or y */
typedef struct {
    Point             *point ; /* x or y                                      */
    VSP_Az               *Az ; /* Ax or Ay                                    */
    VSPINT            weight ; /* c'x or c'y                                  */
} State ;

/* ========================================================================== */
// Constructors and Destructors for Core Types
/* ========================================================================== */

Problem* VSP_construct_problem
(
    VSPINT       num_vertices,
    VSPLONG                nz,
    VSPINT                *Ai, // If NULL, then this is allocated
    VSPLONG               *Ap, // If NULL, then this is allocated
    VSPINT                *Ae, // If NULL, then this is allocated
    VSP_options      *options
) ;

void VSP_destruct_problem(Problem* problem) ;



VSP_Graph *Construct_Graph(VSPINT n, VSPLONG nz, VSPLONG *Ap, VSPINT *Ai, VSPINT*Ae) ;
void Destroy_Graph(VSP_Graph *Graph) ;

VSP_Az* VSP_construct_Az(VSPINT n) ;
void VSP_destruct_Az(VSP_Az* Az) ;

Point* VSP_construct_point(VSPINT n ) ;
void VSP_destruct_point(Point* point) ;

State* VSP_construct_state(VSPINT n);
void VSP_destruct_state(State* state);


/* ========================================================================== */
// Global Variables
/* ========================================================================== */
extern const Problem* G_P ;


#endif
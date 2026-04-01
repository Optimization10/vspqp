#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include<stdbool.h>
#include"VSPQP.h"

typedef struct Stopwatches Stopwatches ;

// #if VSP_INTEGERS_MODE == 0
//     #define VSPINT int ;
//     #define VSPLONG int ;
// #elif VSP_INTEGERS_MODE == 1
//     #define VSPINT int ;
//     #define VSPLONG long long ;
// #elif VSP_INTEGERS_MODE == 2
//     #define VSPINT long long ;
//     #define VSPLONG long long ;
// #else
//     #error "Invalid VSP_INTEGERS_MODE"
// #endif

// typedef double VSPFLOAT ;

typedef enum
{
 VERSION_DEFAULT,
 VERSION_1,
 VERSION_2,
 VERSION_3,
 VERSION_4,
 VERSION_5,
 VERSION_6,
 VERSION_7,
 VERSION_8
} VersionNum ;

typedef enum
{
    STRICT_INCREASE,
    WEAK_INCREASE,
    NO_CHANGE,
    WEAK_DECREASE,
    STRICT_DECREASE,
    UNSPECIFIED_EFFECT
} Obj_Effect ;

typedef enum
{
    N1_OPTIMAL,            // No feasible, strictly beneficial flip
    N1_WEAKLY_OPTIMAL,     // No feasible, strictly beneficial flip AND there is at least one neutral downflip
    N1_STRICTLY_OPTIMAL,   // No feasible, strictly beneficial flip AND there are zero neutral downflips
    NOT_N1_OPTIMAL         // At least one feasible, strictly beneficial flip
} VSP_N1_Optimality ;

/* Usage:

    First, create a VSP_levels_selection_type variable with the desired
    selection type, and then create a call the constructor of VSP_levels
    and store the result in a variable, called, let's say, levels.

    TOP_LEVELS: Then the top (levels->num_top_levels). Note that if this is 0,
                then no level is selected.

    BOT_LEVELS: Then the bottom (levels->num_bot_levels). Note that if this is
                0, then no level is selected.

    ALL_LEVELS: Then all levels are selected.

    NO_LEVELS: Then no levels are selected.

    SPECIFY_SELECTED_LEVELS: Then the levels_selected stack should be filled
                             with the levels to be selected. This stack is
                             added to using the VSP_levels_select function.

    SPECIFY_UNSELECTED_LEVELS: Then the levels_not_selected stack should be
                               filled with the levels to be unselected. This
                               stack is added to using the VSP_levels_unselect
                               function.
*/
typedef enum
{
    TOP_LEVELS,
    BOT_LEVELS,
    ALL_LEVELS,
    NO_LEVELS,
    SPECIFY_SELECTED_LEVELS,
    SPECIFY_UNSELECTED_LEVELS
} VSP_levels_selection_type ;



typedef struct {

    VSPINT num_nodes_total ;
    VSPINT num_edges_total ;

    VSPINT num_nodes_x ;
    VSPINT num_nodes_y ;
    VSPINT num_nodes_sep ;

    VSPFLOAT percent_nodes_x ;
    VSPFLOAT percent_nodes_y ;
    VSPFLOAT percent_nodes_sep ;

    VSPFLOAT time ;
} VSP_solution_stats ;


typedef struct {
    VSPINT  weight_x ;
    VSPINT  weight_y ;
    VSPFLOAT obj_val ;
    VSPINT *x ;
    VSPINT *y ;
} VSP_solve_result ;









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
} VSP_match_type;



/* Store an almost binary point.                                              */
typedef struct {
    VSPINT    *int_part ; /* Integer part*/
    VSPINT     ind_frac ; /* Index of fractional entry. EMPTY if none.        */
    VSPFLOAT   val_frac ; /* Value of fractional entry. No guarantees
                             if ind_frac is EMPTY.                            */
    VSPINT            n ; /* Number of entries in the point                   */
} VSP_point ;

typedef struct
{
    /*--- Generic Info -------------------------------------------------------*/
    VSPINT n  ;      /* Number of vertices */
    VSPLONG nz ;      /* Twice the number of edges */
    char name[200] ; /* Name of the graph */

    /*--- Resolve from Centroid ----------------------------------------------*/
    /* The levels where we tried re-solving */
    VSP_dyn_stk *resolve_tried ;

    /* Stack containing the levels where re-solving a level
       gave a better point than the current best */
    VSP_dyn_stk *resolve_better ;

    /* Stack containing the initial objective at a level */
    VSP_dyn_stk *objIn  ;

    /* Stack containing the final objective at a level */
    VSP_dyn_stk *objOut ;

    /* Stack containing the value of x'Ay */
    VSP_dyn_stk *numPen ;

    /* Stack containing the sum of all the costs in the separator */
    VSP_dyn_stk *sepIn ;

    /* Stack containing the sum of all the costs in the separator */
    VSP_dyn_stk *sepOut ;

    // The number of flips done in the problem
    long long unsigned num_flips ;

} VSP_stats ;



typedef enum
{
    FLIPS_DEC_LEVELS,        /* Specified using a levels struct */
    FLIPS_DEC_NO_COMM,        /* Do flips if the current level was obtained w/o community*/
    FLIPS_DEC_BOUNDARY    /* Do flips if problem->ratioMaxCostNodes > decBoundary */
} VSP_flips_dec_type ;

typedef struct
{
    VSP_flips_dec_type decType ;

    VSP_levels *levels ;
    VSPFLOAT decBoundary ;

} VSP_flips_dec ;

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
    VSP_point         *point ; /* x or y                                      */
    VSP_Az               *Az ; /* Ax or Ay                                    */
    VSPINT            weight ; /* c'x or c'y                                  */
} VSP_state ;


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

typedef enum
{
    UNCOMPRESSED,
    COMPRESSED
} Compression_Status ;

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

/* The algorithm always performs heavy edge matching and brotherly matching.
   This struct records the user's decision about when to enable community
   matching */



void LoggerInit(bool enable);
void LoggerShutdown(void);
void LoggerStart(VSP_state *Xinit, VSP_state *Yinit, VSP_problem *problem, const char *format, ...);
void LoggerEnd(VSP_state *Xfin, VSP_state *Yfin, VSP_problem *problem, const char *format, ...);

#endif
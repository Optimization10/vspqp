#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include "core.h"
#include "VSPQP.h"
#include "time.h"

void VSP_set_centroid
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
) ;

bool Is_Complete_Graph
(
    VSPINT n, VSPLONG nz
) ;

VSPFLOAT VSP_zAw
(
    State *X,
    State *Y
) ;

void Force_To_Binary_Point
(
    State           *Z,
    State           *W,
    Problem   *problem
) ;

VSPINT VSP_penalty_term
(
    VSPINT  *z,
    VSPINT *Aw,
    VSPINT   n
) ;

VSPINT Are_Neighbors
(
    VSPINT vert1,
    VSPINT vert2,
    VSPINT   *Ai,
    VSPLONG *Ap
) ;


void VSP_ResetDumpFile() ;

void VSP_PrintNeighs(VSPINT i, VSPINT *Ai, VSPLONG *Ap) ;

VSPINT VSP_count_nzs(VSPINT *array, VSPINT n) ;

VSPFLOAT VSP_Azi
(
    VSPINT i,
    VSP_Az *Az
) ;

/* Compute c[i] - gamma * Az[i] */
VSPFLOAT VSP_Gradi(VSPINT i, VSPINT *c, VSPFLOAT gamma, VSP_Az *Az) ;


void VSP_GetSepCost
(
    State          *X,
    State          *Y,
    Problem  *problem,
    VSPINT       *sepCost   /* Output*/
) ;

void VSP_copy_array_int
(
    VSPINT *from,
    VSPINT *to,
    VSPINT  n
) ;

VSPINT VSP_bin_points_diff
(
    VSPINT *a,
    VSPINT *b,
    VSPINT n
) ;


void VSP_xAy
(
    State         *X,
    State         *Y,
    Problem *problem,
    VSPFLOAT        *xAy
) ;

VSPLONG VSP_partial_sums_of_array(
    VSPLONG          *arr,
    VSPLONG *partial_sums,
    VSPINT             n
) ;

Problem* VSP_transpose(
    VSPINT           *Ai, /* Row indices of non-zeros  */
    VSPLONG *Ap, /* Column pointers           */
    VSPINT           *Ae, /* Entries                   */
    VSPLONG nz, /* Number of rows/columns    */
    VSPINT             n, /* Number of nonzeros        */
    Options* options
) ;



void VSP_array_to_point(
    VSPFLOAT  *z_float,
    Point       *z
) ;

VSPFLOAT* VSP_Az_struct_to_array (VSP_Az *Az) ;


void VSP_check_against_target_weight (
    VSPINT* x,
    VSPINT ind_frac,
    VSPFLOAT val_frac,
    VSPINT* c,
    VSPINT weight_target,
    VSPINT n
) ;


void Check_Array_Is_Binary
(
    VSPINT *array,
    VSPINT      n
) ;


void VSP_check_linear_constraint (
    VSPINT* int_part,
    VSPINT ind_frac,
    VSPFLOAT val_frac,
    VSPINT* c,
    VSPINT u,
    VSPINT n
) ;

VSPINT VSP_check_problem(
    Problem* problem,
    Options* options
) ;

void Error( 
    char* msg,
    char* file_name,
    VSPINT line_num,
    unsigned long long count
) ;

void Log
(
    const char *filename,
    const char   *format,
    ...
) ;

void VSP_check_weight(State *Z, VSPINT *c) ;

void VSP_check_max_of_array(
    VSPINT        *c,
    VSPINT         n, /* Length of c */
    VSPINT max_given /* Purported maximum value of c */
) ;

void VSP_check_Az(
    Point         *z,
    VSP_Az           *Az,
    Problem *problem
) ;

void VSP_check_arrays_are_equal
(
    const VSPINT *const u, /* Check that vector u equals vector v */
    const VSPINT *const v,
    const VSPINT  n  /* Size of u and v */
) ;

void Check_Arrays_Are_Equal_VSPLONG
(
    const VSPLONG *const u, /* Check that vector u equals vector v */
    const VSPLONG *const v,
    const VSPINT  n  /* Size of u and v */
) ;

void VSP_check_entries_are_between_bounds
(
    VSPINT*  v,
    VSPINT   n, /* Length of v */
    VSPINT  lb,
    VSPINT  up
) ;

void VSP_check_permutation
(
    VSPINT* arr_1,
    VSPINT* arr_2,
    VSPINT n
) ;


void VSP_check_state(State *Z, Problem *problem) ;

void VSP_check_states(State *X, State *Y, Problem *problem) ;

void Check_Part
(
    VSPINT *x,
    VSPINT *Ax,
    VSPINT *y,
    VSPINT *Ay,
    VSPINT wgtx,
    VSPINT wgty,
    Problem *problem
) ;



void VSP_check_point(Point *z) ;


void VSP_light_check_Az
(
    VSP_Az *Az
) ;


void VSP_relaxed_to_vertex_separator(
    State* X,
    State* Y,
    Problem* problem
) ;

void VSP_copy_point(
    Point* from,
    Point* to,
    VSPINT n
) ;


int VSP_comparison_increasing(const void *a, const void *b) ;



int VSP_compare_value_index_floats(const void *a, const void *b) ;

VSPINT VSP_compare_val_ind_structs_inc(const void *a, const void *b) ;


void VSP_extract_matrix_col(
    VSPINT           *Ak, /* Array to store the column */
    VSPINT             k, /* Index of column to extract */
    VSPINT           one, /* Representative of 1's */
    VSPINT           *Ai,
    VSPLONG          *Ap
) ;

VSPFLOAT* Point_to_array(Point *z) ;

void VSP_copy_state
(
    State      *from,
    State        *to
) ;

void VSP_compute_weight(Point *Z, VSPINT *weight_z, VSPINT *c) ;

int VSP_comparison_increasing(const void *a, const void *b) ;

void VSP_resolve_frac_entry(State *Z, State *W, Problem *problem) ;

VSPINT Is_Array_Member
(
    VSPINT candidate,
    VSPINT    *array,
    VSPINT     start,
    VSPINT       end
) ;


bool VSP_equal
(
    VSPFLOAT a,
    VSPFLOAT b
) ;

/* Check that a <= b */
bool VSP_lteq
(
    VSPFLOAT a,
    VSPFLOAT b
) ;



Options* VSP_construct_options();
void VSP_destruct_options(Options* options);

VSPINT** VSP_construct_explicit_matrix (
    VSPINT* values, VSPINT* Ai, VSPLONG nz, VSPLONG* Ap, VSPINT n
) ;

void VSP_destruct_explicit_matrix( VSPINT** A, VSPINT n ) ;
                                        

VSPFLOAT VSP_compute_LP_objective(
    Point    *z,
    VSPINT       *c,
    VSPFLOAT  gamma,
    VSP_Az      *Aw
) ;

double Graph_Density(VSPINT n, VSPLONG nz) ;

void VSP_mat_vec_mult
(
    const VSPINT *const    z,
    VSPINT       *const   Az,
    const VSPINT *const   Ai,
    const VSPLONG *const   Ap,
    const VSPINT           n
) ;

void Compute_Az
(
    Point              *z,
    VSP_Az                *Az,
    VSPINT                *Ai,
    VSPLONG               *Ap,
    VSPINT                  n
) ;

void VSP_Az_pm_col_of_A
(
    VSPINT             *Az, // vector to which we add/subtract
    VSPINT               i, // column to add/subtract
    VSPINT          operation, // 1 if add, -1 if subtract
    Problem   *problem
) ;

VSPFLOAT VSP_compute_objective_0
(
    State *X,
    State *Y,
    Problem *problem,
    VSPFLOAT gamma
) ;

void Compute_ObjValSimple
(
    VSPINT               *x,
    VSPINT              *Ax,
    VSPINT             wgtx,
    VSPINT               *y,
    VSPINT              *Ay,
    VSPINT             wgty,
    Problem    *problem,
    VSPINT       *retObjVal
) ;

VSPFLOAT VSP_compute_obj
(
    State *X,
    State *Y,
    Problem *problem,
    VSPFLOAT gamma
) ;

void Compute_ObjValDbg
(
    VSPINT               *x,
    VSPINT              *Ax,
    VSPINT               *y,
    VSPINT              *Ay,
    Problem    *problem,
    VSPINT       *retObjVal
) ;




void DumpPoint(Point *z, char *outName) ;
void DumpPoints(Point *x, Point *y, char *outName) ;
void DumpIntArray(VSPINT *array, VSPINT len, char *outName) ;
void DumpGrad(VSPINT *c, VSPFLOAT gamma, VSP_Az *Az, char *outName) ;

void VSP_SetFracEntry
(
    VSPINT         indFrac,
    VSPFLOAT       valFrac,
    State           *Z,
    Problem   *problem
) ;

VSPFLOAT VSP_ComputeFracEntry
(
    VSPINT   u,
    VSPINT wgt,
    VSPINT  ci
) ;



void VSP_ComputeGradFlt
(
    VSPFLOAT  *grad,
    VSPINT       *c,
    VSPFLOAT  gamma,
    VSPFLOAT     *z,
    VSPINT      *Ai,
    VSPLONG *Ap,
    VSPINT        n
) ;

VSPFLOAT VSP_LP_solver_0
(
    VSPFLOAT *x,         /* solution                                         */
    VSPFLOAT *d,         /* objective                                        */
    VSPINT *c,           /* constraint                                       */
    VSPINT l,            /* Lower bound. Must be 0 or 1                      */
    VSPINT u,            /* RHS of constraint. Must be at least 1            */
    VSPINT n,            /* size of problem                                  */
    Options *options /* Choose randomness or not                         */
) ;


bool Is_N1_Optimal_Dbg
(
    State         *X,
    State         *Y,
    Problem *problem,
    bool          strict  // If true, check for strict N1 Optimality
) ;

void limited_printf(const char *fmt, ...) ;


bool Is_N1_Optimal_Dbg
(
    State                   *X,
    State                   *Y,
    Problem           *problem,
    bool          strictly_optimal  // If true, check for strict N1 Optimality
) ;


VSPFLOAT VSP_minreduce
(
    State         *X,
    State         *Y,
    Options *options,
    Problem *problem
) ;

void VSP_set_random_initial_guess
(
    State         *X,
    State         *Y,
    Problem *problem,
    Options *options
) ;


VSPINT Count_Ones
(
    VSPINT *x, // assumed binary
    VSPINT  n  // size of x
) ;


void Copy_Data
(
    VSPINT *xfrom,
    VSPINT *Axfrom,
    VSPINT *wgtxfrom,
    VSPINT *yfrom,
    VSPINT *Ayfrom,
    VSPINT *wgtyfrom,
    VSPINT *xto,
    VSPINT *Axto,
    VSPINT *wgtxto,
    VSPINT *yto,
    VSPINT *Ayto,
    VSPINT *wgtyto,
    VSPINT n
) ;


void Get_Permutation(VSPINT *arr, VSPINT n) ;

void Check_Feasibility
(
    VSPINT      *xInt, // integer part of x 
    VSPINT   indFracX, // ignored if EMPTY
    VSPFLOAT valFracX,
    VSPINT      *yInt, // integer part of y
    VSPINT   indFracY, // ignored if EMPTY
    VSPFLOAT valFracY,
    VSPINT         *c,
    VSPINT          u,
    VSPINT          n
) ;



void Check_Az_Int
(
    const VSPINT *const  z,
    const VSPINT *const Az,
    const VSPINT *const Ai,
    const VSPLONG *const Ap,
    const VSPINT         n
) ;

void Check_Weight
(
    const VSPINT    *const z,
    const VSPINT     indFrac, // Ignored if EMPTY
    const VSPFLOAT   valFrac, // Ignored if EMPTY
    const VSPINT    *const c,
    const VSPINT         wgt,
    const VSPINT           n
) ;


void Print_Neighs
(
    VSPINT vtx, VSPINT *Ai, VSPLONG *Ap
) ;



bool Is_Coarsest_Problem(Problem *problem) ;
bool Is_Finest_Problem(Problem *problem) ;
void Swap_Mem(void *a, void *b, size_t size) ;


// Return a pointer to a state containing a deep copy of z, Az, and wgtz
void Save_State_Copy
(
    VSPINT *z, VSPINT *Az, VSPINT wgtz, VSPINT n, State **stateHome
) ;

// Return a pointer to a state containing a deep copy of z, Az, and wgtz
void Save_State_Copy
(
    VSPINT *z, VSPINT *Az, VSPINT wgtz, VSPINT n, State **stateHome
) ;

void Call_Stack_Print() ;


void Fix_Lower_Bound_Issue
(
    State         *Z,
    State         *W,
    Problem *problem
) ;


// base ^ expopnent where both base and exponent are
// 0^0 is 1
VSPINT pow_ints(VSPINT base, VSPINT exponent) ;

bool Implies(bool A, bool B) ;

void IO_Dump_Partition_int32
(
    VSPINT         *x,
    VSPINT         *y,
    VSPINT          n,
    char   *file_name  // Null terminated string. If NULL, then use default name
) ;

bool Is_Singleton
(
    VSPINT    i,
    VSPLONG *Ap
) ;


bool Is_Less_Than(VSPFLOAT a, VSPFLOAT b, VSPFLOAT tol) ;
bool Is_Greater_Than(VSPFLOAT a, VSPFLOAT b, VSPFLOAT tol) ;

void State_Round_Down
(
    State         *Z,
    Problem *problem
) ;

VSPINT Get_Max_Degree(VSPLONG *Ap, VSPINT n) ;

VSPINT Problem_Get_N_Finest(const Problem *P) ;
VSPLONG Problem_Get_NZ_Finest(const Problem *P) ;

// problem.hager6 n nz Ai Ap c x y (integer parts) Ax Ay (integer parts)
// problem.txt Contains u, and any fractional entries of X and Y
void Dump_Problem
(
    const Problem *P,
    const State   *X,
    const State   *Y
) ;

VSPINT State_Get_Wgt_Int_Part
(
    State *Z,
    VSPINT    *c
) ;

VSPINT Get_Rand_Array_Entry
(
    VSPINT *arr,
    VSPINT    n  // size of arr
) ;

#endif
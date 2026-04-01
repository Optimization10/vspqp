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
















/* Store an almost binary point.                                              */
typedef struct {
    VSPINT    *int_part ; /* Integer part*/
    VSPINT     ind_frac ; /* Index of fractional entry. EMPTY if none.        */
    VSPFLOAT   val_frac ; /* Value of fractional entry. No guarantees
                             if ind_frac is EMPTY.                            */
    VSPINT            n ; /* Number of entries in the point                   */
} VSP_point ;





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



typedef enum
{
    UNCOMPRESSED,
    COMPRESSED
} Compression_Status ;


/* The algorithm always performs heavy edge matching and brotherly matching.
   This struct records the user's decision about when to enable community
   matching */



void LoggerInit(bool enable);
void LoggerShutdown(void);
void LoggerStart(VSP_state *Xinit, VSP_state *Yinit, VSP_problem *problem, const char *format, ...);
void LoggerEnd(VSP_state *Xfin, VSP_state *Yfin, VSP_problem *problem, const char *format, ...);

#endif
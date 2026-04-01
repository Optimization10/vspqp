/* Header file to be included only in main.c */

#include "shared_api.h"
#include "other_api.h"
#include "io_api.h" 
#include "VSP_read_problem.h"
#include <time.h>
#include <string.h>
#include <locale.h>
#include "ds_api.h"
#include <stdio.h>





typedef struct
{
    VSPINT numNodesx ;
    VSPINT numNodesy ;
    VSPINT costSep ;
    long double time ; // Time it took

} METIS_Solution_Stats ;

static void SetCustomOptions(VSP_options *options) ;
static void SaveRndSeed(time_t seed) ;
static void PrintSumOfResults
(
    VSPINT           num_sep,
    VSPINT           sum_imb,
    long double         time,
    unsigned int        seed,
    VSPINT        num_graphs
) ;

void VSP_Solve_With_Metis
(
    char                      *name,
    VSPINT                      *Ai,
    VSPLONG *Ap,
    VSPINT                       *c, // if NULL, then it's all 1's
    VSPINT                        n,
    VSPFLOAT                ubRatio,
    long long unsigned         seed,
    METIS_Solution_Stats       *MSS // Ignored if NULL
) ;

void VSP_Print_Performance
(

) ;

VSPINT Prepend
(
    const char      *bufferA,
    char            *bufferB,
    size_t       sizeBufferB
) ;


VSP_solve_result* VSP_solve
(
    VSP_problem *problem,
    VSP_options *options
) ;


void VSP_create_dot_ovi_file
(
    VSP_problem *problem,
    char           *name,    /* From cmnd line. Ex: data / **** / bcspwr08.mtx */
    char         *outDir
) ;


VSP_problem *Read_hager6_File(char *fileName, VSP_options *options) ;


VSP_problem *Read_hager6_File(char *path, VSP_options *options) ;
VSPINT GetNextLine(const char* fileName, char* output, size_t outputSize) ;

#define MAX_HEIGHT 10
#define MAX_WIDTH  255

#define ABS_DIFF(a, b) ( ((a) > (b)) ? ( (a) - (b) ) : ( (b) - (a) ) )

double G_ROLLBACK_TIME = 0 ;


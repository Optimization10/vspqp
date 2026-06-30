#include "VSPQP.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// #include <cstdlib>
#include "time.h"
#include <sys/stat.h>

#define MAX_MATRIX_PATH_LENGTH 1000
#define MAX_LEN_GRAPH_NAME 1000

#define ABS_DIFF(a, b) ( ((a) > (b)) ? ( (a) - (b) ) : ( (b) - (a) ) )
#define CAST(arg, type) ( ((type) (arg)) )
#define TRUE  1
#define FALSE 0

/* The cleanup is done in-place */
void cleanupGraphName(char *graph_name)
{
    /* Remove directories */
    char *last_slash = strrchr(graph_name, '/');
    if (last_slash != NULL)
    {
        memmove(graph_name, last_slash + 1, strlen(last_slash + 1) + 1);
    }

    /* Remove extension */
    char *last_dot = strrchr(graph_name, '.');
    if (last_dot != NULL)
    {
        *last_dot = '\0';
    }

    return ;
}

static void getColTitle
(
    char           *titleGiven,
    char         *titleToWrite
)
{
    if(*titleGiven == '\0')
    {
        strcpy(titleToWrite, "Ovi ") ;
    }
    else
    {
        if( strchr(titleGiven, ',') != NULL)
        {
            printf("\nError: No comma allowed in column title.\n") ;
            exit(1) ;
        }
        strcpy(titleToWrite, titleGiven) ;
    }
    return ;
}


int Directory_Exists(const char *path)
{
    struct stat sb ;
    return stat(path, &sb) == 0 && S_ISDIR(sb.st_mode) ;
}

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
} solution_stats ;



static void Print_Aggregate_Results
(
    VSPINT           num_sep,
    VSPINT           sum_imb,
    long double         time,
    unsigned int        seed,
    VSPINT        num_graphs
)
{
    printf("\nNum Graphs  : %llu" , CAST(num_graphs, long long unsigned)) ;
    printf("\nSep  : %llu" , CAST(num_sep, long long unsigned)) ;
    printf("\nImb  : %llu" , CAST(sum_imb, long long unsigned)) ;
    printf("\nTime : %.2Lf\n\n", time) ;

    FILE *fp = fopen("output/aggregate.txt", "w") ;
    if(fp == NULL)
    {
        printf("\nError: Could not open aggregate file") ;
        return ;
    }
    
    fprintf(fp, "Num Graphs  : %llu" , CAST(num_graphs, long long unsigned)) ;
    fprintf(fp, "\nSep  : %llu" , CAST(num_sep, long long unsigned)) ;
    fprintf(fp, "\nImb  : %llu" , CAST(sum_imb, long long unsigned)) ;
    fprintf(fp, "\nTime : %.2Lf", time) ;
    fprintf(fp, "\nSeed : %llu" , CAST(seed, long long unsigned)) ;
    

    fclose(fp) ;
    return ;
}

solution_stats* compute_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    VSPINT             n,
    VSPLONG           nz,
    VSPFLOAT        time
)
{
    /* Local vars */
    VSPINT i ;
    VSPINT num_nodes_x = 0 ;
    VSPINT num_nodes_y = 0 ;

    /* Allocate memory */
    solution_stats *stats = calloc(1, sizeof(*stats)) ;

    /* Compute number of nodes in x and y */
    for(i = 0 ; i < n ; i++) {
        num_nodes_x += x[i] ;
        num_nodes_y += y[i] ;
    }

    /* Fill out number of nodes in stats struct */
    stats->num_nodes_total = n ;
    stats->num_nodes_x     = num_nodes_x ;
    stats->num_nodes_y     = num_nodes_y ;
    stats->num_nodes_sep   = n - num_nodes_x - num_nodes_y ;

    stats->num_edges_total = nz / 2 ;

    /* Fill out percentages in stats struct */
    stats->percent_nodes_x   = (VSPFLOAT) num_nodes_x * 100 / n ;
    stats->percent_nodes_y   = (VSPFLOAT) num_nodes_y * 100 / n ;
    stats->percent_nodes_sep = (VSPFLOAT) stats->num_nodes_sep * 100 / n ;

    /* Fill out the time */
    stats->time = time ;

    return stats ;
}


void reset_file(char* name)
{
    FILE* file = fopen(name, "w") ;
    fclose(file) ;
    return ;
}


void destruct_solve_result(solve_result *result)
{
    free(result->x) ;
    free(result->y) ;
    free(result) ;

    return ;
}




/******************************************************************************/
// Read Hager6 File
/******************************************************************************/
// Input:  a file to open (absolute or relative path. Prefer absolute)
void Read_hager6_File
(
    char       *path,
    Options *options,

    // Output
    VSPINT    *n,
    VSPINT  **Ai,
    VSPLONG **Ap
)
{

    FILE *fp = fopen(path, "rb") ;
    if (!fp)
    {
        printf("\nERROR: Could not open file: %s", path) ;
        fflush(stdout) ;
        // return NULL ;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Read n and nz
    int32_t ntemp ;
    int64_t nztemp ;
    unsigned long tmp ;
    tmp = fread(&ntemp, sizeof(int32_t), 1, fp) ;
    tmp = fread(&nztemp, sizeof(int64_t), 1, fp) ;
    (void) tmp ; // silence unused result warning

    int32_t *tempAi = malloc(nztemp * sizeof(int32_t)) ;
    if (!tempAi)
    {
        *Ai = NULL ;
        *Ap = NULL ;
        return ;
    } ;

    tmp = fread(tempAi, sizeof(int32_t), nztemp, fp);
    (void) tmp ; // silence unused result warning

    int64_t *tempAp = malloc((ntemp + 1) * sizeof(int64_t));
    if (!tempAp)
    {
        free(tempAi) ;
        *Ai = NULL ;
        *Ap = NULL ;
        return ;
    }  ;

    tmp = fread(tempAp, sizeof(int64_t), ntemp + 1, fp);
    (void) tmp ; // silence unused result warning

    clock_gettime(CLOCK_MONOTONIC, &end) ;
    double read_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9 ;

    fclose(fp);

    // Write n
    *n = ntemp ;

    // Convert Ai array to VSPINTs
    *Ai = malloc(nztemp * sizeof(*Ai));
    if (!Ai)
    {
        free(tempAi);
        free(tempAp);
        // return NULL ;
    }
    for (int64_t i = 0; i < nztemp; i++) (*Ai)[i] = tempAi[i];

    // Convert Ap array to VSPLONGs
    *Ap = malloc((ntemp + 1) * sizeof(*Ap) );
    if (!Ap)
    {
        free(Ai) ;
        free(tempAi) ;
        free(tempAp) ;
        // return NULL ;
    }
    for (int32_t i = 0 ; i <= ntemp ; i++) (*Ap)[i] = tempAp[i] ;

    // Free memory
    free(tempAi);
    free(tempAp);
}



// Input: A file where each line has at most 1000 chars

// Output: on the nth call, it returns the nth line
// On any calls after the last line has been read, the function returns NULL
// Does NOT prepend a path
// Writes the next line from the list file to the provided buffer
// If it ever reads an empty line, 
VSPINT GetNextLine(const char* fileName, char* output, size_t outputSize)
{
    // Static variables to maintain state between calls
    static FILE* fp = NULL; // File pointer for reading the file
    static VSPINT eof_reached = 0; // Flag to track EOF state

    // If EOF was reached previously, return 0
    if (eof_reached) return 0;

    // Open file if not already open
    if (fp == NULL)
    {
        // Open the file
        fp = fopen(fileName, "r");
        if (fp == NULL) return 0;
    }

    // Temporary buffer for reading a line
    char line[1001];
    // Read the next line from the file
    if (fgets(line, sizeof(line), fp) != NULL)
    {
        // Remove trailing newline
        line[strcspn(line, "\n")] = '\0';

        // Copy line to output buffer
        strncpy(output, line, outputSize - 1);

        // Ensure null-termination
        output[outputSize - 1] = '\0';

        // Return 1 to indicate success
        return 1;
    }
    else
    {
        // Close file and set EOF state
        fclose(fp);
        fp = NULL;
        eof_reached = 1;
        // Return 0 to indicate end of file
        return 0;
    }
}



// Prepend the null terminated string in buffer A to the null terminated
// string in buffer B
// Returns 1 if success, and 0 if failure (if the combined length of the strings
// exceeds sizeBufferB)
VSPINT Prepend
(
    const char      *bufferA,
    char            *bufferB,
    size_t       sizeBufferB
)
{
    size_t lenA = strlen(bufferA) ;
    size_t lenB = strlen(bufferB) ;

    if (lenA + lenB + 1 > sizeBufferB)
    {
        return 1 ;
    }

    // Shift bufferB to the right by lenA positions
    memmove(bufferB + lenA, bufferB, lenB + 1) ;

    // Copy bufferA to the beginning
    memcpy(bufferB, bufferA, lenA) ;

    return 1 ; // success
}



/******************************************************************************/

void Write_Results
(
    char           *graph_name,
    solution_stats   *solStats,
    Options           *options
)
{
    static long long unsigned count = 0 ;
    count++ ;

    char cpy_graph_name[MAX_LEN_GRAPH_NAME] ;
    strcpy(cpy_graph_name, graph_name) ;
    cleanupGraphName(cpy_graph_name) ;

    FILE *file ;
    if(count == 1)
    {
        // Reset file and write column titles
        file = fopen("output/results.csv", "w") ;
        fprintf(file, "id,n,edges,sep,imb,time\n" ) ;
    }
    else
    {
        file = fopen("output/results.csv", "a") ;
    }

    /* Write data */
    fprintf(file, "%s,", cpy_graph_name) ;
    fprintf(file, "%lld,", (long long) solStats->num_nodes_total) ;
    fprintf(file, "%lld,", (long long) solStats->num_edges_total) ;
    fprintf(file, "%lld,", (long long) solStats->num_nodes_sep) ;
    fprintf(file, "%lld,", (long long) labs((long long) solStats->num_nodes_x - (long long) solStats->num_nodes_y)) ;
    fprintf(file, "%.6f\n", solStats->time) ;

    fclose(file) ;
    return ;
}


int main(int argc, char** argv)
{
    /*------------------------------------------------------------------------*/
    // User Controls for Choosing a Problem Set
    /*------------------------------------------------------------------------*/
     
    // Text file from which to read matrix names, one per line
    char matricesList[] = "data/to-solve" ;

    // Set data directory
    char matricesDir[256];
    if (Directory_Exists("/home/ovi/data/suitesparse/symmetric/") )
    {
        strcpy(matricesDir, "/home/ovi/data/suitesparse/symmetric/") ;
    }
    else
    {
        strcpy(matricesDir, "/media/hager/datafile/symmetric/") ;
    }

    // strcpy(matricesDir, "/home/onechita/network-repository/") ;
    // strcpy(matricesDir, "/home/ovi/data/network-repository/networks/") ;

    char *name ;
    
    
    /*------------------------------------------------------------------------*/
    // Reset progress file and error log
    reset_file("./output/progress.txt") ;
    reset_file("./output/errors.txt") ;
    unsigned int seed ;

    // seed = (unsigned VSPINT) time(NULL) ;
    seed = 1 ;

    srand(seed) ;
    printf("\nSeed : %u\n\n", seed ) ;

    /* Local vars */
    VSPINT                 sum_sep ;
    VSPINT                 sum_imb ;
    long double           sum_time ;
    struct timespec     start, end ;
    long double       elapsed_time ;
    Options           *options ;
    solve_result       *result = NULL ;


    /*** Options **************************************************************/
    options = construct_options() ;
    set_default_options(options) ;

    // Custom options
    // options->sort_cols     = true ;
    options->attempt_second_start = FALSE ;
    // options->u_ratio_coarsest  = 0.45 ;
    options->u_ratio_coarsest  = 0.48 ;
    // options->u_ratio_coarsest  = 0.42 ;
    options->initial_guess = bisection_initial_guess ;
    // options->ub_based_on_bisection = true ;
    options->initial_project = FALSE ;
    // options->initial_project = TRUE ;
    options->FM_Limit = 320 ;
    options->try_both_covers = FALSE ;
    options->bisection_ngp = 3 ;
    options->vsp_grad_proj = FALSE ;
    options->vsp_ngp = 3 ;

    Print_Options(options, "output/opts.txt") ;
    /*** I/O ******************************************************************/
    /* 1. Set column title for results.csv */
    char filePath[MAX_MATRIX_PATH_LENGTH] = "" ;

    VSPINT num_graphs = 0 ;
    double total_time = 0 ;

    /*--------------------------------------------------------------------*/
    // Solve all problems using flips
    /*--------------------------------------------------------------------*/
    sum_imb  = 0 ;
    sum_sep  = 0 ;
    sum_time = 0 ;

    VSPINT n, *Ai ;
    VSPLONG *Ap ;
    while(GetNextLine(matricesList, filePath, MAX_MATRIX_PATH_LENGTH))
    {
        num_graphs++ ;
        srand(seed) ;

        // Read the problem
        Prepend(matricesDir, filePath, MAX_MATRIX_PATH_LENGTH) ;
        Read_hager6_File(filePath, options, &n, &Ai, &Ap) ;

        // Log the current problem name
        name = strrchr(filePath, '/') + 1 ;
        FILE *f = fopen("./output/progress.txt", "a") ;
        fprintf(f, "Starting problem: %s\n", name) ;
        fclose(f) ;

        // Solve the problem
        restart_decide_reset(options->restart_decision) ;
        clock_gettime(CLOCK_REALTIME, &start) ;
        result = VSP_Solve(n, Ap, Ai, options, name) ;
        clock_gettime(CLOCK_REALTIME, &end);


        // Compute elapsed and aggregate time
        elapsed_time = compute_elapsed_time(start, end) ;
        total_time += elapsed_time ;


        // Compute solution stats
        solution_stats *sltn_stats =
        compute_solution_stats(result->x, result->y, n, Ap[n], elapsed_time);
        sum_sep += sltn_stats->num_nodes_sep ;
        sum_imb += ABS_DIFF(sltn_stats->num_nodes_x, sltn_stats->num_nodes_y) ; 
        sum_time += elapsed_time ;

        // Write to output/results.csv
        Write_Results(filePath, sltn_stats, options) ;

        // Free memory
        free(Ai) ;
        free(Ap) ;
        destruct_solve_result(result) ;
        free(sltn_stats) ;
    }
    Print_Aggregate_Results(sum_sep, sum_imb, sum_time, seed, num_graphs) ;


    /* Free memory */
    destruct_options(options) ;
    return 0 ;
}




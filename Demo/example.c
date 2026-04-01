#include "VSPQP.h"
// #include "main.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// #include <cstdlib>
#include "time.h"
#include <sys/stat.h>

#define MAX_MATRIX_PATH_LENGTH 1000

struct timespec  G_START_TIME, G_END_TIME ;
long double G_COMPRESS_TIME = 0, G_FLIPS_TIME = 0, G_COARSEN_TIME = 0,
G_HEM_TIME = 0, G_BM_TIME = 0 ;

int Directory_Exists(const char *path)
{
    struct stat sb ;
    return stat(path, &sb) == 0 && S_ISDIR(sb.st_mode) ;
}


VSP_solution_stats* VSP_compute_solution_stats
(
    VSPINT            *x,
    VSPINT            *y,
    VSP_problem *problem,
    VSPFLOAT        time
)
{
    /* Local vars */
    VSPINT i ;
    VSPINT n = problem->Graph->n ;
    VSPINT num_nodes_x = 0 ;
    VSPINT num_nodes_y = 0 ;

    /* Allocate memory */
    VSP_solution_stats *stats = calloc(1, sizeof(*stats)) ;

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

    stats->num_edges_total = problem->Graph->nz / 2 ;

    /* Fill out percentages in stats struct */
    stats->percent_nodes_x   = (VSPFLOAT) num_nodes_x * 100 / n ;
    stats->percent_nodes_y   = (VSPFLOAT) num_nodes_y * 100 / n ;
    stats->percent_nodes_sep = (VSPFLOAT) stats->num_nodes_sep * 100 / n ;

    /* Fill out the time */
    stats->time = time ;

    return stats ;
}


void VSP_reset_file(char* name)
{
    FILE* file = fopen(name, "w") ;
    fclose(file) ;
    return ;
}

static void SaveRndSeed(time_t seed)
{
    FILE* file = fopen("output/VSP_error_log.txt", "w") ;
    fprintf(file, "Seed: %li", seed) ;
    fclose(file) ;

    file = fopen("output/random_seed.txt", "w") ;
    fprintf(file, "Seed: %li", seed) ;
    fclose(file) ;
    return ;
}

void VSP_reset_output_files()
{
    VSP_reset_file("output/VSP_results.txt") ;
    VSP_reset_file("output/VSP_error_log.txt") ;
    VSP_reset_file("output/VSP_stats.txt") ;
    return ;
}


/******************************************************************************/
// Read Hager6 File
/******************************************************************************/
// Input:  a file to open (absolute or relative path. Prefer absolute)
// Output: problem struct. Obtained by calling constructor, freeing its Ai,
//         Ap arrays, and setting them to the Ai and Ap arrays set in this
//         program
VSP_problem *Read_hager6_File(char *path, VSP_options *options)
{

    FILE *fp = fopen(path, "rb") ;
    if (!fp)
    {
        printf("\nERROR: Could not open file: %s", path) ;
        fflush(stdout) ;
        return NULL ;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Read n and nz
    int32_t ntemp ;
    int64_t nztemp ;
    fread(&ntemp, sizeof(int32_t), 1, fp) ;
    fread(&nztemp, sizeof(int64_t), 1, fp) ;

    int32_t *tempAi = malloc(nztemp * sizeof(int32_t)) ;
    if (!tempAi) return NULL ;

    fread(tempAi, sizeof(int32_t), nztemp, fp);

    int64_t *tempAp = malloc((ntemp + 1) * sizeof(int64_t));
    if (!tempAp)
    {
        free(tempAi) ;
        return NULL ;
    }

    fread(tempAp, sizeof(int64_t), ntemp + 1, fp);

    clock_gettime(CLOCK_MONOTONIC, &end) ;
    double read_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9 ;

    fclose(fp);

    // Write the processing time
    double fileSize ;
    struct stat st  ;
    char *fileName = strrchr(path, '/') + 1 ;
    if (stat(path, &st) == 0)
    {
        // Convert bytes to GB
        fileSize = (double)st.st_size / (1024 * 1024 * 1024); 

        fp = fopen("./output/FileReadingTime.csv", "a") ;
        fprintf
        (
            fp, "\n%s,%f,%f,%f",
            fileName, fileSize, read_time,fileSize/read_time
        ) ;
        fclose(fp) ;
    }
    else
    {
        printf("\nCould not read file size") ;
    }

    // Convert Ai array to VSPINTs
    VSPINT *Ai = malloc(nztemp * sizeof(*Ai));
    if (!Ai)
    {
        free(tempAi);
        free(tempAp);
        return NULL ;
    }
    for (int64_t i = 0; i < nztemp; i++) Ai[i] = tempAi[i];

    // Convert Ap array to VSPLONGs
    VSPLONG *Ap = malloc((ntemp + 1) * sizeof(*Ap) );
    if (!Ap)
    {
        free(Ai) ;
        free(tempAi) ;
        free(tempAp) ;
        return NULL ;
    }
    for (int32_t i = 0 ; i <= ntemp ; i++) Ap[i] = tempAp[i] ;

    // Free memory
    free(tempAi);
    free(tempAp);


    // Set up problem to return
    VSPINT  n  = (VSPINT) ntemp ;
    VSPLONG nz = (VSPLONG) nztemp ;


    VSP_problem *problem = VSP_construct_problem(n, Ap[n], options) ;
    free(problem->Graph->Ai) ;
    free(problem->Graph->Ap) ;
    problem->Graph->Ai = Ai ;
    problem->Graph->Ap = Ap ;

    problem->Graph->Ai = Ai ;
    problem->Graph->Ap = Ap ;

    // Write file name
    strcpy(problem->name, fileName) ;

    return problem ;
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

static void SetCustomOptions(VSP_options *options)
{
    /*------------------------------------------------------------------------*/
    /* Coarsening Options                                                     */
    /*------------------------------------------------------------------------*/
    options->use_coarsening             = true ;
    options->max_size_of_coarsest_graph = 10  ;

    options->compress_threshold = 0.7 ;

    /*------------------------------------------------------------------------*/
    // Tim Davis Coarsening with Boundary for Community Matching
    /*------------------------------------------------------------------------*/
    options->matching_strategy          = HEMSR ;
    options->do_community_match         = true ;

    options->coarsenDec = VSP_construct_coarsen_dec(CRSN_DEC_useCommWithBoundary) ;
    options->coarsenDec->decBounday = 0.6 ;

    /*------------------------------------------------------------------------*/
    // Low gamma matching
    /*------------------------------------------------------------------------*/
    // options->matching_strategy = LOW_GAMMA_MATCHING ;



    /*------------------------------------------------------------------------*/
    /* Restart Level Options                                                  */
    /*------------------------------------------------------------------------*/
    options->restart_decision = VSP_construct_restart_decision() ;

    // /* All levels in [2, infty) */
    options->restart_decision->use_levels = true ;
    options->restart_decision->levels     = VSP_construct_levels(SPECIFY_UNSELECTED_LEVELS) ;
    VSP_levels_unselect(0, options->restart_decision->levels) ;
    VSP_levels_unselect(1, options->restart_decision->levels) ;

    /* Using Rule */
    // WARNING!!!! RULE UNSUPPORTED RIGHT NOW
    // options->restart_decision->use_rule = true ;
    // options->restart_decision->restarts_allowed      = 1 ;
    // options->restart_decision->rule_param            = 0.04 ;
    // options->restart_decision->finest_level_allowed  = 2 ;


    // /* Disable restarts */
    // options->restart_decision->use_levels = true ;
    // options->restart_decision->levels = VSP_construct_levels(NO_LEVELS) ;

    /*------------------------------------------------------------------------*/
    /* Gamma Reduction Options                                                */
    /*------------------------------------------------------------------------*/
    options->use_gamma_reduction      = true ;

    options->max_num_vert_small_graph = 8000 ;
    options->num_gr_for_small_graphs  = 4 ;
    options->num_gr_for_large_graphs  = 1 ;


    // trial
    // options->max_num_vert_small_graph = 16000 ;
    // options->num_gr_for_small_graphs  = 2 ;

    // options->max_num_vert_small_graph = 20000 ;
    // options->num_gr_for_small_graphs  = 2 ;



    options->num_its_for_small_graphs = 1 ;
    options->num_its_for_large_graphs = 1 ;


    /*------------------------------------------------------------------------*/
    /* Relax Upper Bound Options                                              */
    /*------------------------------------------------------------------------*/
    options->relaxUB      = true ;
    options->ratioUBLimit = 0.47 ;


    // options->relaxUB      = false ; 
    // options->ratioUBLimit = 0.6 ;

    /*------------------------------------------------------------------------*/
    /* Starting Guess Options                                                 */
    /*------------------------------------------------------------------------*/
    options->initial_guess     = VSP_centroid_initial_guess ;
    
    // best so far
    options->centroidVal       = 0.35 ;

    /*------------------------------------------------------------------------*/
    /* FM Options                                                             */
    /*------------------------------------------------------------------------*/
    // Best so far
    options->FM_Limit = 300 ;

    // options->FM_Limit = 200 ;
    // options->FM_Limit = 250 ;
    // options->FM_Limit = 300 ;
    // options->FM_Limit = 350 ;

    /*------------------------------------------------------------------------*/
    /* Miscellaneous Options                                                  */
    /*------------------------------------------------------------------------*/
    options->sep_at_each_level = false ;


    /* ====================================================================== */
    // Mountain Climbing Options
    /* ====================================================================== */
    options->tol_LP_opt_checker = 1E-10 ;

    /* ====================================================================== */
    // Simple Coarsening options
    /* ====================================================================== */
    options->max_density     = 0.8 ;
    options->limit_num_nodes = 100 ;
    options->compression_threshold = 0.51 ;

    return ;
}


int main(int argc, char** argv)
{
    /*------------------------------------------------------------------------*/
    // User Controls for Choosing a Problem Set
    /*------------------------------------------------------------------------*/
     
    // Text file from which to read matrix names, one per line
    char matricesList[] = "doc/mats/ManualMatsList.txt" ;
    // char matricesList[] = "doc/mats/1.txt" ;
    // char matricesList[] = "doc/mats/int32-n-sorted.txt" ;

    // Set data directory
    char matricesDir[256];
    if (Directory_Exists("/home/ovi/data/symmetric/"))
    {
        strcpy(matricesDir, "/home/ovi/data/symmetric/") ;
    }
    else
    {
        strcpy(matricesDir, "/media/hager/datafile/symmetric/") ;
    }

    // The maximum path length of a problem
    
    /*------------------------------------------------------------------------*/


    // SuiteSparse_start() ;
    VSP_reset_output_files() ; 
    unsigned int seed ;

    // seed = (unsigned VSPINT) time(NULL) ;
    seed = 1 ;

    srand(seed) ;
    SaveRndSeed(seed) ;
    printf("\nSeed : %u", seed ) ;

    /* Local vars */
    VSPINT                 sum_sep ;
    VSPINT                 sum_imb ;
    long double           sum_time ;
    struct timespec     start, end ;
    long double       elapsed_time ;
    VSP_options           *options ;
    VSP_problem           *problem ;
    VSP_solve_result       *result = NULL ;

    /*** Options **************************************************************/
    options = VSP_construct_options() ;
    VSP_set_default_options(options) ;
    SetCustomOptions(options) ;

    /*** I/O ******************************************************************/
    /* 1. Set column title for Numerical_Results.csv */
    char col_title[] = "Sep" ;
    char filePath[MAX_MATRIX_PATH_LENGTH] = "" ;
  

    VSP_reset_file("./output/SolveLog.txt") ;

    VSP_reset_file("./output/FileReadingTime.csv") ;
    FILE *fp = fopen("./output/FileReadingTime.csv", "a") ; 
    fprintf(fp, "File,Size(GB),Time(s),Speed(GB/s)") ;
    fclose(fp) ;


    VSPINT num_graphs = 0 ;
    double total_time = 0 ;


    /*--------------------------------------------------------------------*/
    // Solve all problems using flips
    /*--------------------------------------------------------------------*/
    sum_imb  = 0 ;
    sum_sep  = 0 ;
    sum_time = 0 ;

    while(GetNextLine(matricesList, filePath, MAX_MATRIX_PATH_LENGTH))
    {
        num_graphs++ ;
        srand(seed) ;

        // Append the directory path
        Prepend(matricesDir, filePath, MAX_MATRIX_PATH_LENGTH) ;

        VSP_restart_decide_reset(options->restart_decision) ;

        // /* 2. Read the problem and set default data */
        problem = Read_hager6_File(filePath, options) ;
        VSP_set_default_problem_data(problem) ;

        /* 3. Solve the problem */ 
        FILE *f = fopen("./output/SolveLog.txt", "a") ;
        fprintf(f, "Starting problem: %s\n", problem->name) ;
        fclose(f) ;
        clock_gettime(CLOCK_REALTIME, &start);
        // result = VSP_solve(problem, options) ;
        result = Solve_v2(problem, options) ;
        clock_gettime(CLOCK_REALTIME, &end);

        /* Compute total time for the VSP solver */
        elapsed_time = VSP_compute_elapsed_time(start, end) ;
        total_time += elapsed_time ;


        // /* 4. Compute solution stats */
        VSP_solution_stats *sltn_stats =
        VSP_compute_solution_stats(result->x, result->y, problem, elapsed_time);
        sum_sep += sltn_stats->num_nodes_sep ;
        sum_imb += ABS_DIFF(sltn_stats->num_nodes_x, sltn_stats->num_nodes_y) ; 
        sum_time += elapsed_time ;

        // /* 5. Write results to a csv file */
        VSP_write_to_csv_2(filePath, col_title, sltn_stats, problem, options) ;
        Write_Performance_Profile
        (
            problem,
            sltn_stats->num_nodes_sep,
            sltn_stats->num_nodes_x,
            sltn_stats->num_nodes_y,
            sltn_stats->time
        ) ;

        // /* 6. Free memory */
        VSP_destruct_problem(problem) ;
        VSP_destruct_solve_result(result) ;
        free(sltn_stats) ;
    }
    PrintSumOfResults(sum_sep, sum_imb, sum_time, seed, num_graphs) ;
    VSP_Print_Performance() ;

    /* Free memory */
    VSP_destruct_options(options) ;

    printf("Flips time: %Lf\n", G_FLIPS_TIME) ;
    printf("HEM time: %Lf\n", G_HEM_TIME) ;
    printf("BM time: %Lf\n", G_BM_TIME) ;
    printf("Coarsen time: %Lf\n", G_COARSEN_TIME) ;

    printf("Matching time: %Lf\n", G_Performance_Struct.Matching_Seconds) ;

    VSP_Print_Performance() ;

    return 0 ;
}



static void PrintSumOfResults
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
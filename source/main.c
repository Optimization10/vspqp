#include "main.h"
#include "shared_types.h"
// #include <cstdlib>
#include "time.h"


typedef enum
{
    SOLVE_WITH_METIS,
    SOLVE_WITH_VSPQP
} SolveAlg ;

struct timespec  G_START_TIME, G_END_TIME ;
long double G_COMPRESS_TIME = 0, G_FLIPS_TIME = 0, G_COARSEN_TIME = 0,
G_HEM_TIME = 0, G_BM_TIME = 0 ;

int main(int argc, char** argv)
{
    /*------------------------------------------------------------------------*/
    // User Controls for Choosing a Problem Set
    /*------------------------------------------------------------------------*/
     
    // Text file from which to read matrix names, one per line
    char matricesList[] = "ManualMatsList.txt" ;

    // Absolute or relative path of directory containing the matrices
    // Make sure to include a trailing slash

    // Matrices on ovi computer
    char matricesDir[] = "/home/ovi/data/symmetric/" ;

    // Matrices dir on hager6
    // char matricesDir[] = "/media/hager/datafile/symmetric/" ;

    // Choose between VSPQP and METIS
    SolveAlg solve = SOLVE_WITH_VSPQP ;
    // SolveAlg solve = SOLVE_WITH_METIS ;

    // The maximum path length of a problem
    #define MAX_MATRIX_PATH_LENGTH 1000
    /*------------------------------------------------------------------------*/


    SuiteSparse_start() ;
    VSP_reset_output_files() ;
    unsigned int seed ;

    // seed = (unsigned int) time(NULL) ;
    // seed = 1760026178 ;
    seed = 1760033187 ;
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
    if(solve == SOLVE_WITH_VSPQP)
    {
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
            WriteGlobalMatrixName(filePath) ;

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

            // /* 4. Compute solution stats */
            VSP_solution_stats *sltn_stats =
            VSP_compute_solution_stats(result->x, result->y, problem, elapsed_time);
            sum_sep += sltn_stats->num_nodes_sep ;
            sum_imb += ABS_DIFF(sltn_stats->num_nodes_x, sltn_stats->num_nodes_y) ; 
            sum_time += elapsed_time ;

            // /* 5. Write results to a csv file */
            VSP_write_to_csv_2(filePath, col_title, sltn_stats, problem, options) ;

            // /* 6. Free memory */
            VSP_destruct_problem(problem) ;
            VSP_destruct_solve_result(result) ;
            free(sltn_stats) ;
        }
        PrintSumOfResults(sum_sep, sum_imb, sum_time, seed, num_graphs) ;
        VSP_Print_Performance() ;
    }
    else if(solve == SOLVE_WITH_METIS)
    {
        /*--------------------------------------------------------------------*/
        // Solve all problems with METIS
        /*--------------------------------------------------------------------*/    
        METIS_Solution_Stats MSS ;    
        VSPINT sumSepMetis  = 0 ;
        long double sumTimeMetis = 0 ;
        VSPINT sumImbMetis  = 0 ;
        while( GetNextLine(matricesList, filePath, MAX_MATRIX_PATH_LENGTH) )
        {
            num_graphs++ ;
            srand(seed) ;
            WriteGlobalMatrixName(filePath) ;

            // Append the directory path
            Prepend(matricesDir, filePath, MAX_MATRIX_PATH_LENGTH) ;

            problem = Read_hager6_File(filePath, options) ;

            // Record the problem being solved
            FILE *f = fopen("./output/SolveLog.txt", "a") ;
            fprintf(f, "Starting problem: %s\n", problem->name) ;
            fclose(f) ;

            VSP_set_default_problem_data(problem) ;
            VSP_Solve_With_Metis 
            (
                problem->name, problem->Graph->Ai, problem->Graph->Ap, problem->Graph->c,
                problem->Graph->n , 0.6, seed, &MSS
            ) ;
            sumSepMetis  += MSS.costSep ;
            sumImbMetis  += ABS_DIFF( MSS.numNodesx, MSS.numNodesy) ; 
            sumTimeMetis += MSS.time ;

            VSP_destruct_problem(problem) ;
        }
        printf("\n\nMETIS:") ;
        PrintSumOfResults(sumSepMetis, sumImbMetis, sumTimeMetis, seed, num_graphs);
        /*--------------------------------------------------------------------*/
    }
    else printf("\n Error: Need to select solve alg") ;

    /* Free memory */
    VSP_destruct_options(options) ;

    printf("Flips time: %Lf\n", G_FLIPS_TIME) ;
    printf("HEM time: %Lf\n", G_HEM_TIME) ;
    printf("BM time: %Lf\n", G_BM_TIME) ;
    printf("Coarsen time: %Lf\n", G_COARSEN_TIME) ;

    return 0 ;
}




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

    options->max_num_vert_small_graph = 80 ;
    options->num_gr_for_small_graphs  = 4 ;
    options->num_gr_for_large_graphs  = 1 ;

    options->num_its_for_small_graphs = 1 ;
    options->num_its_for_large_graphs = 1 ;


    /*------------------------------------------------------------------------*/
    /* Relax Upper Bound Options                                              */
    /*------------------------------------------------------------------------*/
    // options->relaxUB      = true ;
    // options->ratioUBLimit = 0.47 ;

    // options->relaxUB      = true ;
    // options->ratioUBLimit = 0.4 ;


    options->relaxUB      = false ; 
    options->ratioUBLimit = 0.6 ;

    /*------------------------------------------------------------------------*/
    /* Starting Guess Options                                                 */
    /*------------------------------------------------------------------------*/
    options->initial_guess     = VSP_centroid_initial_guess ;
    options->centroidVal       = 0.4 ;
    // options->centroidVal       = 0.3 ;

    /*------------------------------------------------------------------------*/
    /* FM Options                                                             */
    /*------------------------------------------------------------------------*/
    options->FM_Limit = 200 ;

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

    return ;
}




static void WriteGlobalMatrixName(char* name)
{
    strcpy(global_matrix, name) ;
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
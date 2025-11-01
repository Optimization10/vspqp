#include "shared_api.h"
#include "metis.h"

typedef struct
{
    VSPINT numNodesx ;
    VSPINT numNodesy ;
    VSPINT costSep ;
    long double time ; // Time it took

} METIS_Solution_Stats ;


//      ASSUMPTION ON UPPER BOUND
//
// In my formulation of the problem, I have a u which represents the maximum
// cost of a partition.
// If I want to communicate this u to METIS, I assume that I have to do
// options[ufactor] = [2 * (u / sum of all ci in graph) - 1] * 1000
// See sfm.c:56 in the METIS library, my notes in Obsidian, and the METIS docs
// I am PRETTY sure this assumption is correct
void VSP_Solve_With_Metis
(
    char                      *name,
    VSPINT                      *Ai,
    VSPLONG *Ap,
    VSPINT                       *c, // if NULL, then it's all 1's
    VSPINT                        n,
    VSPFLOAT                ubRatio,
    VSPINT                     seed,
    METIS_Solution_Stats       *MSS // Ignored if NULL
)
{
    static long long unsigned count = 0 ;
    count++ ;

    /* Convert from data from VSPINT to idx_t which is used by metis */
    idx_t numVertices   = (idx_t) n ;
    idx_t numNonzeros   = (idx_t) (Ap[n]) ;
    idx_t *adjPointers  = malloc((numVertices + 1) * sizeof(*adjPointers)) ;
    idx_t *adjIndices   = malloc( numNonzeros * sizeof(*adjIndices) ) ;
    idx_t *costVertices = malloc(numVertices * sizeof(*costVertices)) ;

    idx_t i, sumAllci = 0 ;
    for(i = 0 ; i < numVertices + 1 ; i++) adjPointers [i] = (idx_t) Ap[i] ;
    for(i = 0 ; i < numNonzeros     ; i++) adjIndices  [i] = (idx_t) Ai[i] ;
    for(i = 0 ; i < numVertices     ; i++)
    {
        costVertices[i] = (c == NULL) ? 1 : (idx_t) c [i] ;
        sumAllci       += costVertices[i] ;
    }


    /* Allocate memory for partition result */
    idx_t *partition = malloc(numVertices * sizeof(*partition)) ;

    /* Set metis options */
    idx_t options[METIS_NOPTIONS] ;
    METIS_SetDefaultOptions(options) ;

    // Coarsening Scheme
    // options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM ;
    // options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM ;

    // Initial Partitioning Scheme
    options[METIS_OPTION_IPTYPE] = 
    // METIS_IPTYPE_GROW;
    // METIS_IPTYPE_RANDOM;
    // METIS_IPTYPE_EDGE;
    METIS_IPTYPE_NODE;
    // METIS_IPTYPE_METISR8;
    
    // Refinement Scheme
    // options[METIS_OPTION_RTYPE] = 
    // METIS_RTYPE_FM;       
    // METIS_RTYPE_GREEDY;    
    // METIS_RTYPE_SEP2SIDED; 
    // METIS_RTYPE_SEP1SIDED; 
    
    /* Unbalance Factor */    
    idx_t ufactor =  199 ;
    // idx_t ufactor = (idx_t) ( (2.0 * ubRatio - 1) * 1000 ) ;
    options[METIS_OPTION_UFACTOR] =  ufactor ;

    /* Random Seed */    
    options[METIS_OPTION_SEED] = (idx_t) seed ;
    // options[METIS_OPTION_SEED] = 2 ;
    // options[METIS_OPTION_SEED] = rand() % 100000 ;

    /* Advanced coarsening techniques */
    //Aggregate nodes which have same adjacency lists.
    // options[METIS_OPTION_COMPRESS] = 0 ;
    //Match nodes using 2-hops when standard matching fails
    // options[METIS_OPTION_NO2HOP] = 1 ;
    // options[METIS_OPTION_NO2HOP] = 0 ;

    /* Number of refinement iterations. Extremely weird, seems to make almost no difference */
    options[METIS_OPTION_NITER] = 10 ;
    // options[METIS_OPTION_NITER] = 1000000 ;

    /* Number of separators computed by METIS at each level */
    // options[METIS_OPTION_NSEPS] = 1 ;
    // options[METIS_OPTION_NSEPS] = 2 ;
    // options[METIS_OPTION_NSEPS] = 1 ;


    // Compute vertex separator
    idx_t sepSize ;
    struct timespec start, end ;
    clock_gettime(CLOCK_REALTIME, &start);
    UNUSED VSPINT metisExit = METIS_ComputeVertexSeparator
                        (
                            &numVertices, adjPointers, adjIndices,
                            costVertices, options, &sepSize, partition
                        ) ;
    clock_gettime(CLOCK_REALTIME, &end);
    long double elapsed_time = VSP_compute_elapsed_time(start, end) ;


    // Compute x, y vectors and the cost of each set
    VSPINT *x   = calloc(numVertices, sizeof(*x)) ;
    VSPINT *y   = calloc(numVertices, sizeof(*y)) ;
    VSPINT *Ay  = calloc(numVertices, sizeof(*Ay)) ;
    VSPINT *sep = calloc(numVertices, sizeof(*sep)) ;


    VSPINT costx = 0, costy = 0, costsep = 0, numvertx = 0, numverty = 0,
           numvertsep = 0 ;

    for(i = 0 ; i < numVertices ; i++)
    {
        if(partition[i] == 0)
        {
            x[i] = 1 ;
            numvertx++ ;
            costx += c[i] ;
        }
        else if(partition[i] == 1)
        {
            y[i] = 1 ; 
            numverty++ ;
            costy += c[i] ;
        }
        else if(partition[i] == 2)
        {
            numvertsep++ ;
            costsep += c[i] ;
        }
        else
        {
            printf("\n\nInvalid entry in METIS partition array\n\n") ;
        }
    }

    // Check that the two separator sizes match
    if(sepSize != costsep)
    {
        printf("\n\nSep sizes don't match in METIS\n\n") ;
    }

    // Check that there are no edges between x and y
    VSP_mat_vec_mult(y, Ay, Ai, Ap, n) ;
    for(i = 0 ; i < numVertices ; i++)
    {
        if( x[i] && (y[i] || Ay[i]) )
        {
            printf("\n\nNot a valid vertex separator\n\n") ;
            break ;
        }
    }


    FILE *file ;
    if(count == 1)
    {
        file = fopen("output/METIS_Numerical_Results.csv", "w") ;

        /* Write column titles */
        fprintf(file, "Graph, ") ;
        fprintf(file, "Vert, ") ;
        fprintf(file, "Edges, ") ;
        fprintf(file, "Total_Graph_Cost, ") ;
        fprintf(file, "METIS_Sep, ") ;
        fprintf(file, "METIS_Time, ") ;
        fprintf(file, "Num_Vert_x, ") ;
        fprintf(file, "Cost_x, ") ;
        fprintf(file, "Rel_Cost_x, ") ;
        fprintf(file, "Num_Vert_y, ") ;
        fprintf(file, "Cost_y, ") ;
        fprintf(file, "Rel_Cost_y, ") ;
        fprintf(file, "Cost_Sep, ") ;
        fprintf(file, "Rel_Cost_Sep, ") ;
        fprintf(file, "Upper_Bal, ") ;
        fprintf(file, "Rel_Upper_Bal\n") ;

    }
    else
    {
        file = fopen("output/METIS_Numerical_Results.csv", "a") ;
    }


    // Calculate relative costs, ensuring no division by zero
    double rel_cost_x = (sumAllci != 0)   ? ((double)costx / sumAllci) : 0.0;
    double rel_cost_y = (sumAllci != 0)   ? ((double)costy / sumAllci) : 0.0;
    double rel_cost_sep = (sumAllci != 0) ? ((double)costsep / sumAllci) : 0.0;
    double rel_upper_bal = ubRatio;

    // Write a row to the csv file
    fprintf(file, "%s, ", name);
    fprintf(file, "%llu, ", (unsigned long long)numVertices);
    fprintf(file, "%llu, ", (unsigned long long)numNonzeros);
    fprintf(file, "%llu, ", (unsigned long long)sumAllci);
    fprintf(file, "%llu, ", (unsigned long long)sepSize);
    fprintf(file, "%.6Lf, ", elapsed_time);
    fprintf(file, "%llu, ", (unsigned long long)numvertx);
    fprintf(file, "%llu, ", (unsigned long long)costx);
    fprintf(file, "%.6f, ", rel_cost_x);
    fprintf(file, "%llu, ", (unsigned long long)numverty);
    fprintf(file, "%llu, ", (unsigned long long)costy);
    fprintf(file, "%.6f, ", rel_cost_y);
    fprintf(file, "%llu, ", (unsigned long long)costsep);
    fprintf(file, "%.6f, ", rel_cost_sep);
    fprintf(file, "%llu, ", (unsigned long long)ufactor);
    fprintf(file, "%.6f\n", rel_upper_bal); // End the line with newline

    fclose(file) ;

    // ------------------- Fill Out the METIS_Solution_Stats Struct -------------------

    // Assign costs
    if(MSS != NULL)
    {
        MSS->costSep = costsep;
        MSS->time    = elapsed_time ;
        MSS->numNodesx = costx ;
        MSS->numNodesy = costy ;
    }
    // --------------------------------------------------------------------------------



    // Free memory
    free(adjPointers) ;
    free(adjIndices) ;
    free(costVertices) ;
    free(partition) ;
    free(sep) ;
    free(x) ;
    free(y) ;
    free(Ay) ;

    return ;
}


// METIS options from James. The ones with only one set of slashes are
// the ones selected by James
    // // Coarsening Scheme
    // options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM ;
    // // options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM ;

    // // Initial Partitioning Scheme
    // options[METIS_OPTION_IPTYPE] = 
    // //METIS_IPTYPE_GROW;
    // //METIS_IPTYPE_RANDOM;
    // //METIS_IPTYPE_EDGE:
    // METIS_IPTYPE_NODE;
    // //METIS_IPTYPE_METISR8;
    
    // // Refinement Scheme
    // // options[METIS_OPTION_RTYPE] = 
    // // METIS_RTYPE_FM;       
    // // METIS_RTYPE_GREEDY;    
    // METIS_RTYPE_SEP2SIDED; 
    // // METIS_RTYPE_SEP1SIDED; 
    
    // /* Unbalance Factor */    
    // idx_t ufactor =  199 ;
    // // idx_t ufactor = (idx_t) ( (2.0 * ubRatio - 1) * 1000 ) ;
    // options[METIS_OPTION_UFACTOR] =  ufactor ;

    // /* Random Seed */    
    // options[METIS_OPTION_SEED] = 1 ;   
    // // options[METIS_OPTION_SEED] = 2 ;   
    // // options[METIS_OPTION_SEED] = seed;   

    // /* Advanced coarsening techniques */
    // //Aggregate nodes which have same adjacency lists.
    // options[METIS_OPTION_COMPRESS] = 0 ;
    // //Match nodes using 2-hops when standard matching fails
    // options[METIS_OPTION_NO2HOP] = 1 ;

    // /* Number of refinement iterations */
    // options[METIS_OPTION_NITER] = 10 ;
    // // options[METIS_OPTION_NITER] = 100 ;

    // /* Number of separators computed by METIS at each level */
    // options[METIS_OPTION_NSEPS] = 1 ;
    // // options[METIS_OPTION_NSEPS] = 40 ;

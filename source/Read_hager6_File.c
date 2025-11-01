#include "shared_api.h"
#include "string.h"
#include <inttypes.h>
#include <stdio.h>
#include <sys/stat.h>



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
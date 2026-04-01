#include "core.h"
#include "logger.h"
#include "shared_api.h"
#include <stdarg.h>
#include "dynamic_stack.h"

struct Logger_Record
{
    VSP_state         *X ;
    VSP_state         *Y ;
    VSP_problem *problem ;
    char            *msg ;
} ;

// Function to create formatted string buffer
static char* LoggerCreateFormattedBuffer(const char *format, va_list args)
{
    if (!format)
    {
        char *buffer = (char *)malloc(1);
        buffer[0] = '\0'; // Empty string for no format
        return buffer;
    }

    va_list args_copy;
    va_copy(args_copy, args);
    
    // Calculate required buffer size
    VSPINT size = vsnprintf(NULL, 0, format, args_copy) + 1;
    va_end(args_copy);
    
    char *buffer = (char *)malloc(size);
    if (!buffer) {
        return NULL; // Handle allocation failure
    }
    
    vsnprintf(buffer, size, format, args);
    return buffer;
}

void LoggerInit(bool enable)
{
    LOG_G.enabled = enable ;

    if(!enable) return ;
    LOG_G.log = fopen("output/log.txt", "w") ;
    LOG_G.initRecords = DynamicStack_Init() ;
}

void LoggerDestroy()
{
    if(!LOG_G.enabled) return ;
    DynamicStack_Destroy(LOG_G.initRecords) ;
    fclose(LOG_G.log) ;
}

void LoggerStart
(
    VSP_state *X, VSP_state *Y, VSP_problem *problem,
    const char *format, ...
)
{
    if(!LOG_G.enabled) return ;
    static long long unsigned counter = 0 ;



    // Handle optional format string
    va_list args;
    va_start(args, format);
    char *msg = LoggerCreateFormattedBuffer(format, args);
    va_end(args);

    Logger_Record *rec = calloc(1, sizeof(*rec)) ;
    VSP_state *Xcpy = VSP_construct_state(problem->Graph->n) ;
    VSP_state *Ycpy = VSP_construct_state(problem->Graph->n) ;
    VSP_copy_state(X, Xcpy) ;
    VSP_copy_state(Y, Ycpy) ;
    rec->X = Xcpy ;
    rec->Y = Ycpy ;
    rec->problem = problem ;
    rec->msg = msg ;

    VSPINT size = DynamicStack_GetSize(LOG_G.initRecords) ;
    fprintf(LOG_G.log, "\n") ;
    for(VSPINT i = 0 ; i < size ; i++)
    {
        fprintf(LOG_G.log, "    ") ;
    }
    fprintf(LOG_G.log, "Start : %s", msg) ;

    DynamicStack_Push(LOG_G.initRecords, (void*) rec) ;
    return ;
}

void LoggerStop(VSP_state *X, VSP_state *Y, VSP_problem *problem)
{
    if(!LOG_G.enabled) return ;
    Logger_Record *rec = (Logger_Record*) DynamicStack_Pop(LOG_G.initRecords) ;    
    
    VSPFLOAT objInit, objFin, imbInit, imbFin, gamma ;
    gamma = CAST(problem->Graph->max_vtx_cost, VSPFLOAT) ;
    objInit = VSP_compute_obj(rec->X, rec->Y, rec->problem, gamma) ;
    objFin  = VSP_compute_obj(X, Y, problem, gamma) ;

    VSPINT size = DynamicStack_GetSize(LOG_G.initRecords) ;

    fprintf(LOG_G.log, "\n") ;
    for(VSPINT i = 0 ; i < size ; i++)
    {
        fprintf(LOG_G.log, "    ") ;
    }
    fprintf(LOG_G.log, "Stop  : %s   Obj change %%: %f", rec->msg, 100*(objFin - objInit)/objFin ) ;

    VSP_destruct_state(rec->X) ;
    VSP_destruct_state(rec->Y) ;
    free(rec->msg) ;
    free(rec) ;
    return ;
}
#ifndef SHARED_GLOBALS_H
#define SHARED_GLOBALS_H


#include "shared_types.h"

extern const VSP_problem* G_P ;

/******************************************************************************
    GLOBAL PERFORMANCE MEASURING STRUCT
*******************************************************************************/



    #define TIMER_INIT \
    struct timespec hidden_timespec_struct_start ; \
    struct timespec hidden_timespec_struct_end ;

    #define TIMER_START \
    struct timespec hidden_timespec_struct_start ; \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_start)

    #define TIMER_ELAPSED \
    (VSP_compute_elapsed_time(hidden_timespec_struct_start, hidden_timespec_struct_end))
    
    #define TIMER_END \
    struct timespec hidden_timespec_struct_end ;   \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_end)

    #define PERFORMANCE_UPDATE_TIME(PREFIX) \
    // do { \
    //     G_Performance_Struct.PREFIX##_Count++; \
    //     G_Performance_Struct.PREFIX##_Seconds += \
    //     (   (hidden_timespec_struct_end.tv_sec - hidden_timespec_struct_start.tv_sec) \
    //     + (hidden_timespec_struct_end.tv_nsec - hidden_timespec_struct_start.tv_nsec) / 1e9L ) ;\
    // } while(0)


#endif
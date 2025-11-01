#ifndef SHARED_GLOBALS_H
#define SHARED_GLOBALS_H

#include "shared_types.h"
extern const VSP_problem* G_P ;

/******************************************************************************
    GLOBAL PERFORMANCE MEASURING STRUCT
*******************************************************************************/
#if MEASURE_PERFORMANCE
    typedef struct
    {
        long double All_Coarsening_Seconds ;
        long long unsigned All_Coarsening_Count ;

        long double MC_Seconds ;
        long long unsigned MC_Count ;

        long double FM_Seconds ;
        long long unsigned FM_Count ;

        long double Flips_Seconds ;
        long long unsigned Flips_Count ;

        long double FM_Flips_Seconds ;
        long long unsigned FM_Flips_Count ;

        long double Relaxed_To_VSP_Seconds ;
        long long unsigned Relaxed_To_VSP_Count ;

        long double Matching_Seconds ;
        long long unsigned Matching_Count ;
    } VSP_Performance ;

    extern VSP_Performance G_Performance_Struct ;

    #define TIMER_START \
    struct timespec hidden_timespec_struct_start ; \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_start)

    #define TIMER_END \
    struct timespec hidden_timespec_struct_end ;   \
    clock_gettime(CLOCK_REALTIME, &hidden_timespec_struct_end)

    #define PERFORMANCE_UPDATE_TIME(PREFIX) \
    do { \
        G_Performance_Struct.PREFIX##_Count++; \
        G_Performance_Struct.PREFIX##_Seconds += \
        (   (hidden_timespec_struct_end.tv_sec - hidden_timespec_struct_start.tv_sec) \
        + (hidden_timespec_struct_end.tv_nsec - hidden_timespec_struct_start.tv_nsec) / 1e9L ) ;\
    } while(0)
#else
    #define TIMER_START
    #define TIMER_END
    #define PERFORMANCE_UPDATE_TIME(PREFIX)
#endif

#endif
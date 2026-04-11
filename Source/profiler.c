/**
 * @file profiler.c
 * @brief Profile code and write results to CSV
 *
 * Accumulated timings are written as CSV rows via Timer_Flush()
 * and reset after each flush. Use the TIMER_START / TIMER_STOP macros
 * defined in profiler.h to instrument individual phases.
 */

#include <stdio.h>
#include "profiler.h"
#include "core.h"

/** Human-readable names for each timer phase, indexed by timer_id. */
static const char *timer_names[] = {"flips", "mc", "fm", "coarsen", "total", "equiv"};

/** Cumulative wall-clock seconds for each timer phase. */
double timings[T_COUNT];

/** Number of times each timer phase has been stopped. */
int counts[T_COUNT];

/**
 * Zero all accumulated timings and call counts.
 */
static void Timer_Reset(void)
{
    for (int i = 0; i < T_COUNT; i++)
    {
        timings[i] = 0.0;
        counts[i] = 0;
    }
}

/**
 * Write one CSV row per active timer phase, then reset all timers.
 *
 * Phases with zero calls are silently skipped. If no CSV file is open
 * the call is a no-op.
 *
 * @param problem  Problem name written to the CSV (NULL becomes "NA").
 * @param level    Coarsening level written to the CSV.
 */
void Timer_Flush(VSPINT level)
{
    FILE* fp ;
    static int first_call = 1 ;
    if(first_call)
    {
        first_call = 0 ;
        fp = fopen("output/timing.csv", "w") ;
        fprintf(fp, "problem,level,phase,total_s,calls,avg_s\n");
    }
    else fp = fopen("output/timing.csv", "a") ;

    for (int i = 0; i < T_COUNT; i++)
    {
        if (counts[i] > 0)
            fprintf(fp, "%s,%d,%s,%.6f,%d,%.6f\n",
                    G_Current_Problem_Name, level, timer_names[i],
                    timings[i], counts[i], timings[i] / counts[i]);
    }
    Timer_Reset();
    fclose(fp) ;
}
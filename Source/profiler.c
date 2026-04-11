/**
 * @file profiler.c
 * @brief Profile code and write results to CSV
 *
 * Accumulated timings are written as CSV rows via VSP_timer_flush()
 * and reset after each flush. Use the TIMER_START / TIMER_STOP macros
 * defined in profiler.h to instrument individual phases.
 */

#include <stdio.h>
#include "profiler.h"

/** Human-readable names for each timer phase, indexed by timer_id. */
static const char *timer_names[] = {"partition", "coarsen", "refine", "total"};

/** Output CSV file handle, opened by VSP_timer_open(). */
static FILE *csv;

/** Cumulative wall-clock seconds for each timer phase. */
double timings[T_COUNT];

/** Number of times each timer phase has been stopped. */
int counts[T_COUNT];

/**
 * Open (or reopen) the profiling CSV file and write the header row.
 *
 * If a file is already open it is closed first. On failure the file
 * handle is left NULL and subsequent flushes become no-ops.
 *
 * @param filename  Path to the output CSV file.
 */
void VSP_timer_open(const char *filename)
{
    if (csv) fclose(csv);
    csv = fopen(filename, "w");
    if (!csv) return ;
    fprintf(csv, "problem,level,phase,total_s,calls,avg_s\n");
}

/**
 * Zero all accumulated timings and call counts.
 */
void VSP_timer_reset(void)
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
void VSP_timer_flush(const char *problem, int level)
{
    if (!csv) return ;
    for (int i = 0; i < T_COUNT; i++)
    {
        if (counts[i] > 0)
            fprintf(csv, "%s,%d,%s,%.6f,%d,%.6f\n",
                    problem ? problem : "NA", level, timer_names[i],
                    timings[i], counts[i], timings[i] / counts[i]);
    }
    VSP_timer_reset();
}

/**
 * Close the CSV file, if open.
 */
void VSP_timer_close(void)
{
    if (csv)
    {
        fclose(csv) ;
        csv = NULL ;
    }
}
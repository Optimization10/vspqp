#include <stdio.h>
#include "profiler.h"

static const char *timer_names[] = {"partition", "coarsen", "refine", "total"};
static FILE *csv;

double timings[T_COUNT];
int counts[T_COUNT];

void VSP_timer_open(const char *filename) {
    if (csv) fclose(csv);
    csv = fopen(filename, "w");
    if (!csv) { perror("VSP_timer_open"); return; }
    fprintf(csv, "problem,level,phase,total_s,calls,avg_s\n");
}

void VSP_timer_reset(void) {
    for (int i = 0; i < T_COUNT; i++) {
        timings[i] = 0.0;
        counts[i] = 0;
    }
}

void VSP_timer_flush(const char *problem, int level) {
    if (!csv) return;
    for (int i = 0; i < T_COUNT; i++) {
        if (counts[i] > 0)
            fprintf(csv, "%s,%d,%s,%.6f,%d,%.6f\n",
                    problem ? problem : "NA", level, timer_names[i],
                    timings[i], counts[i], timings[i] / counts[i]);
    }
    VSP_timer_reset();
}

void VSP_timer_close(void) {
    if (csv) { fclose(csv); csv = NULL; }
}


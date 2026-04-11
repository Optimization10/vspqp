#ifndef PROFILER_H
#define PROFILER_H

#include <time.h>

enum timer_id { T_PARTITION, T_COARSEN, T_REFINE, T_TOTAL, T_COUNT };

extern double timings[T_COUNT];
extern int counts[T_COUNT];

#define TIMER_START(id) \
    struct timespec _ts_##id; \
    clock_gettime(CLOCK_MONOTONIC, &_ts_##id)

#define TIMER_STOP(id) do { \
    struct timespec _te_##id; \
    clock_gettime(CLOCK_MONOTONIC, &_te_##id); \
    timings[id] += (_te_##id.tv_sec - _ts_##id.tv_sec) \
                 + (_te_##id.tv_nsec - _ts_##id.tv_nsec) * 1e-9; \
    counts[id]++; \
} while(0)

void VSP_timer_open(const char *filename);
void VSP_timer_flush(const char *problem, int level);
void VSP_timer_close(void);
void VSP_timer_reset(void);

#endif
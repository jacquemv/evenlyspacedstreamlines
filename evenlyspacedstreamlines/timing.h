#include <omp.h>
#include <stdio.h>

static double start_times[32];
static double cumulated_times[32] = {0, 0, 0, 0,  0, 0, 0, 0, 
                                     0, 0, 0, 0,  0, 0, 0, 0,
                                     0, 0, 0, 0,  0, 0, 0, 0, 
                                     0, 0, 0, 0,  0, 0, 0, 0};

void timing_on(int n)
{
    start_times[n] = omp_get_wtime();
}

void timing_off(int n)
{
    cumulated_times[n] += omp_get_wtime() - start_times[n] ;
}

void timing_print(int n, const char* name)
{
    fprintf(stderr, "%-20s [%.3f s]\n", name, cumulated_times[n]);
}

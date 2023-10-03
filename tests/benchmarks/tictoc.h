#include <omp.h>

double start_time = 0;

void tic()
{
    start_time = omp_get_wtime();
}

void toc(const char* name)
{
    printf("%30s [%.3f s]\n", name, omp_get_wtime() - start_time);
    tic();
}

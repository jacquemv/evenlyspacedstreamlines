
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "tictoc.h"
#include "maxsubinterval.cpp"


void random_data(int n, int* lo, int* hi, int k)
{
    for (int i=0;i<n;i++) {
        if (rand() % k > 0) {
            lo[i] = 0;
            hi[i] = n-1;
        } else {
            if (i == 0) lo[i] = 0; else lo[i] = rand() % i;
            if (i == n-1) hi[i] = n-1; else hi[i] = i + rand() % (n-i);
        }
    }
}

void bench(int order, int n, int iter_max, int k)
{
    int *lo = new int[n];
    int *hi = new int[n];
    int iter = 0;
    int a, b;
    while (iter++ < iter_max) {
        if (k < 0) {
            for (int i=0;i<n;i++) {
                lo[i] = 0;
                hi[i] = n-1;
            }
        } else
            random_data(n, lo, hi, k);
        int a, b;
        switch (order) {
            case 3:
                max_admissible_subinterval3(n, lo, hi, a, b);
                break;
            case 2:
                max_admissible_subinterval2(n, lo, hi, a, b);
                break;
            case 1:
                max_admissible_subinterval1(n, lo, hi, a, b);
                break;
            case 0:
                max_admissible_subinterval(n, lo, hi, a, b, 120000);
        }
    }
    delete [] lo;
    delete [] hi;
}

int main()
{
    int n = 10000, niter = 10000;
    tic();
    bench(2, n, niter, 50);
    toc("meth2");
    bench(1, n, niter, 50);
    toc("meth1");
    bench(0, n, niter, 50);
    toc("meth0");
}
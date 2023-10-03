
#include <stdio.h>
#include <stdlib.h>
#include <cassert>

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

void check_consistency(int n, int iter_max, int k)
{
    int *lo = new int[n];
    int *hi = new int[n];
    int iter = 0;

    while (iter++ < iter_max) {
        random_data(n, lo, hi, k);
        int a1, b1, a2, b2, a3, b3;
        max_admissible_subinterval3(n, lo, hi, a3, b3);
        max_admissible_subinterval2(n, lo, hi, a2, b2);
        max_admissible_subinterval1(n, lo, hi, a1, b1);
        assert(a1 == a2 && b1 == b2);
        assert(a3 == a2 && b3 == b2);
        assert(is_subinterval_admissible(n, lo, hi, a1, b1));
    }

    delete [] lo;
    delete [] hi;
}

int main()
{
    check_consistency(1, 5, 1);
    check_consistency(2, 5, 1);
    check_consistency(3, 5, 1);
    check_consistency(20, 1000, 1);
    check_consistency(100, 1000, 1);
    check_consistency(100, 1000, 3);
    check_consistency(100, 1000, 10);
}
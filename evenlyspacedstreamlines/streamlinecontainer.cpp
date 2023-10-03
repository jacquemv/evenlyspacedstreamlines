#include <stdio.h>
#include "algebra.h"
#include "geometry.h"
#include "streamlinecontainer.h"

//-----------------------------------------------------------------------------
StreamlineContainer::StreamlineContainer()
{
    nmax = 0;
    points = NULL;
    idx = NULL;
    n = 0;

}

//-----------------------------------------------------------------------------
StreamlineContainer::StreamlineContainer(int nb_points)
{
    nmax = nb_points;
    assert( nmax >= 2 );
    points = new double [3*nmax];
    idx = new int [nmax-1];
    n = 0;
}
//-----------------------------------------------------------------------------
StreamlineContainer::~StreamlineContainer()
{
    delete [] points;
    delete [] idx;
}

//-----------------------------------------------------------------------------
void StreamlineContainer::push_segment(double* P, double* Q, int face)
{
    if (n == 0) {
        vcopy(points, P);
        vcopy(points+3, Q);
        idx[0] = face;
        n = 2;
    } else if (n == 2) {
        double *first = points, *second = points+3;
        int code = connect_segments(first, second, P, Q);
        if (code % 2 == 0) vcopy(points+6, Q); else vcopy(points+6, P);
        if (code / 2 == 0) vswap(first, second);
        idx[1] = face;
        n = 3;
    } else {
        assert( n < nmax && n > 0);
        double* last = points+3*n-3;
        double* next = last + 3;
        if (vdist2(last, P) < vdist2(last, Q))
            vcopy(next, Q);
        else
            vcopy(next, P);
        idx[n-1] = face;
        n++;
    }
}

//-----------------------------------------------------------------------------
int StreamlineContainer::first_vector(double* vector)
{
    vdiff(vector, points, points+3);
    return idx[0];
}

//-----------------------------------------------------------------------------
void StreamlineContainer::reverse()
{
    #define SWAP(type, x, y) { type tmp; tmp = y; y = x; x = tmp; }
    for (int i=0;i<(n-1)/2;i++) {
        int j = n-2 - i;
        SWAP(int, idx[i], idx[j-1]);
    }
    for (int i=0;i<n/2;i++) {
        int j = n-1 - i;
        SWAP(double, points[3*i], points[3*j]);
        SWAP(double, points[3*i+1], points[3*j+1]);
        SWAP(double, points[3*i+2], points[3*j+2]);
    }
    #undef SWAP
}

//-----------------------------------------------------------------------------
double StreamlineContainer::length()
{
    double L = 0;
    for (int i=0;i<n-1;i++) {
        L += vdist(points+3*i, points+3*i+3);
    }
    return L;
}

//-----------------------------------------------------------------------------
void StreamlineContainer::print()
{
    for (int i=0;i<n;i++) {
        printf("%.5f %.5f %.5f  ", points[3*i], points[3*i+1], points[3*i+2]);
    }
}
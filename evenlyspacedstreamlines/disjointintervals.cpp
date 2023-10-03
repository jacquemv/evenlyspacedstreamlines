#include "disjointintervals.h"

//-----------------------------------------------------------------------------
DisjointIntervals::DisjointIntervals()
{
	nmax = 4;
	n = 0;
	cov = 0.0;
	t = new double [2*nmax];
}

//-----------------------------------------------------------------------------
DisjointIntervals::DisjointIntervals(int default_nmax)
{
	nmax = default_nmax;
	n = 0;
	cov = 0.0;
	t = new double [2*nmax];
}

//-----------------------------------------------------------------------------
DisjointIntervals::~DisjointIntervals()
{
	delete [] t;
}

//-----------------------------------------------------------------------------
void DisjointIntervals::reallocate(int new_nmax)
{
    nmax = new_nmax;
    double* tnew = new double [2*nmax];
    for (int i=0;i<2*n;i++) tnew[i] = t[i];
    delete [] t;
    t = tnew;
}

//-----------------------------------------------------------------------------
inline void DisjointIntervals::clear()
{
	n = 0;
	cov = 0.0;
}

//-----------------------------------------------------------------------------
inline bool DisjointIntervals::is_empty()
{
	return n == 0;
}

//-----------------------------------------------------------------------------
inline bool DisjointIntervals::is_full()
{
	return cov >= 1;
}

//-----------------------------------------------------------------------------
inline bool DisjointIntervals::contains(double x)
{
	// branchless version
	bool b = 0;
	for (int i=0;i<n;i++) b ^= (x >= t[2*i]) ^ (x > t[2*i+1]);
	return b;
	/* equivalent version
	for (int i=0;i<n;i++) {
		if ((x >= t[2*i]) && (x <= t[2*i+1])) return 1;
	}
	return 0;*/
}

//-----------------------------------------------------------------------------
inline void DisjointIntervals::update_coverage()
{
	cov = 0;
	for (int i=0;i<n;i++) {
		cov += t[2*i+1]-t[2*i];
	}
}

//-----------------------------------------------------------------------------
inline double DisjointIntervals::coverage()
{
	return cov;
}

//-----------------------------------------------------------------------------
void DisjointIntervals::insert(double a, double b)
{
	if (b < a) return;
	if (a < 0) a = 0;
	if (b > 1) b = 1;
	// insert first interval
	if (n == 0) {
		n = 1;
		t[0] = a;
		t[1] = b;
		cov = b-a;
		return;
	}
    // find position of a and b
	int ia = 0;
    while ((ia < 2*n) && (a > t[ia])) ia++;
    int ib = ia;
    while ((ib < 2*n) && (b >= t[ib])) ib++;
    // [a, b] lies in a free interval
	if (ia == ib) {
		if (ia % 2 == 1) return;
		n++;
		// reallocate memory
		if (n > nmax) reallocate(2*nmax);
		// insert interval
		for (int i=2*n-1;i>ia+1;i--) t[i] = t[i-2];
		t[ia] = a;
		t[ia+1] = b;
		cov += b-a;
		return;
	}
	// extend interval
	int ma = ia / 2;
	if (ia % 2 == 0)
		t[2*ma] = a;
	if (ib % 2 == 0)
		t[2*ma+1] = b;
	else
		t[2*ma+1] = t[ib];
	// remove intervals
	int mb = (ib-1) / 2;
	int d = mb - ma;
	if (d > 0) {
		for (int i=2*ma+2;i<2*n-2*d;i++) t[i] = t[i+2*d];
		n -= d;
	}
	update_coverage();
}

//-----------------------------------------------------------------------------
double DisjointIntervals::pick_random(int rnd)
{
	if (is_full()) return -1.0;
	double x = ((double)rnd/(double)(RAND_MAX));
	if (is_empty())
		return x;
	x *= 1-cov;
	if (x < t[0]) return x;
	x -= t[0];
	for (int i=1;i<n;i++) {
		double L = t[2*i] - t[2*i-1];
		if (x < L)
			return x + t[2*i-1];
		x -= L;
	}
	return x + t[2*n-1];	
}

//-----------------------------------------------------------------------------
double DisjointIntervals::pick_random()
{
	unsigned int rnd;
	if (is_full()) return -1.0;
	if (is_empty()) {
		//#pragma omp critical
		rnd = rand();
		return (double)rnd/(double)(RAND_MAX);
	}
	//#pragma omp critical
	rnd = rand();
	double x = ((double)rnd/(double)(RAND_MAX)) * (1-cov);
	if (x < t[0]) return x;
	x -= t[0];
	for (int i=1;i<n;i++) {
		double L = t[2*i] - t[2*i-1];
		if (x < L)
			return x + t[2*i-1];
		x -= L;
	}
	return x + t[2*n-1];	
}
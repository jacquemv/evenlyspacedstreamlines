#ifndef DISJOINTINTERVALS_H_
#define DISJOINTINTERVALS_H_

class DisjointIntervals {
public:
    // set of disjoint closed intervals
    // S = [t(0), t(1)[ u [t(2), t(3)] u ... u [t(2n-2), t(2n-1)]
    // where 0 <= t(0) < t(1) < ... < t(2n-1) <= 1
	int n;
	double* t; // bounds of intervals (size 2n)
	
	DisjointIntervals(); // default: allocate for 4 intervals
	DisjointIntervals(int default_nmax);
	~DisjointIntervals(); 

	inline void clear(); // set n = 0
	inline bool is_empty(); // test if the set is empty
    inline bool is_full(); // test if the set = [0, 1]
	inline double coverage(); // returns the coverage (between 0 and 1)

    inline bool contains(double x); // test is x is in one of the intervals

    void insert(double a, double b); // add an interval [a, b] and reduce to a 
                                     // set of disjoint intervals

	double pick_random(); // select a random number in [0, 1] but ouside the 
                          // intervals; returns -1.0 if coverage is complete
	double pick_random(int rnd); // same but provide a precomputed
								 // random number

private:
    int nmax; // nmax = maximum number of intervals allocated
	double cov; // coverage = t(1)-t(0) + ... + t(2n-2)-t(2n-1)

	void reallocate(int new_nmax);
	inline void update_coverage(); // internal use
};

#endif
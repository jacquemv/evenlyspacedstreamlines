/*  For each i, 0 <= i < n, a range [lo[i], hi[i]] is defined such that
        0 <= lo[i] <= i <= hi[i] < n
    A subinterval [a, b] included in [0, n-1] is said to be admissible if
        lo[i] <= a <= b <= hi[i] for all i in [a, b]
    The objective is to find the longest admissible subinterval.

    Interpretation: i is a segment of a polygonal line; 
    [lo[i], hi[i]] is a consecutive set of segments that satisfy some
    constraint with respect to i; we are looking for a consecutive set of
    segments that satisfy all constraints

    Algorithms 1, 2 and 3 have a worst-case complexity of O(n), O(n^2) and 
    O(n^3) respectively. 
    Algorithm 2 is generally faster unless n is very large.
*/

#include <deque>
#include <cassert>

//-----------------------------------------------------------------------------
bool is_subinterval_admissible(int n, int* lo, int* hi, int a, int b)
{
    for (int k=a;k<=b;k++) {
        if (lo[k] > a || hi[k] < b) return 0;
    }
    return 1;
}

//-----------------------------------------------------------------------------
void max_admissible_subinterval3(int n, int* lo, int* hi, int &a, int &b)
{
    a = b = 0;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        if (is_subinterval_admissible(n, lo, hi, i, j) && j-i > b-a) {
            a = i;
            b = j;
        }
    }
}

//-----------------------------------------------------------------------------
void max_admissible_subinterval2(int n, int* lo, int* hi, int &a, int &b)
{
    a = b = 0;
    int i = 0, j = 0;
    int lo_max, hi_min;

    for (i=0;i<n;i++) {
        lo_max = lo[i], hi_min = hi[i];
        j = i;
        while (j < hi_min) {
            if (lo[j+1] > i) break;
            j++;
            if (hi[j] < hi_min) hi_min = hi[j];
            if (lo[j] > lo_max) lo_max = lo[j];
        }
        if (j-i > b-a) {
            a = i;
            b = j;
        }
    }
}

//-----------------------------------------------------------------------------
class SlidingWindowMinMax{
public:
	std::deque<int> deque_min;
	std::deque<int> deque_max;
    int left, right;
    int *xmin, *xmax;

    SlidingWindowMinMax(int* values_for_min, int* values_for_max) {
        xmin = values_for_min;
        xmax = values_for_max;
		deque_min.push_back(xmin[0]);
		deque_max.push_back(xmax[0]);
        left = 0;
        right = 0;
    }

	inline int min() {
		return deque_min.front();
	}

	inline int max() {
		return deque_max.front();
	}

    bool check_max() {
        int m = xmax[left];
        for (int i=left+1;i<=right;i++) if (xmax[i] > m) m = xmax[i];
        return m == max();
    }

    bool check_min() {
        int m = xmin[left];
        for (int i=left+1;i<=right;i++) if (xmin[i] < m) m = xmin[i];
        return m == min();
    }

    inline bool test() {
        return (max() <= left) && (min() >= right);
    }

    inline bool test_push() {
        return (max() <= left) && (min() >= right+1) && 
               (xmax[right+1] <= left) && (xmin[right+1] >= right+1);
    }

    inline int length() {
        return right-left+1;
    }

	void push_right() {
        right++;

        int val = xmax[right];
		while (!deque_max.empty() && val > deque_max.back())
			deque_max.pop_back();
		deque_max.push_back(val);

        val = xmin[right];
		while (!deque_min.empty() && val < deque_min.back())
			deque_min.pop_back();
		deque_min.push_back(val);
	}
	
	void pop_left() {
        int val = xmax[left];
		assert(val <= deque_max.front());
        if (val == deque_max.front()) deque_max.pop_front();

        val = xmin[left];
		assert(val >= deque_min.front());
        if (val == deque_min.front()) deque_min.pop_front();

        left++;
	}
};

//-----------------------------------------------------------------------------
void max_admissible_subinterval1(int n, int* lo, int* hi, int &a, int &b)
{
    a = b = 0;
    if (n <= 1) return;
    SlidingWindowMinMax minmax(hi, lo);
    while (minmax.left < n) {
        while (minmax.right+1 < n && minmax.test_push()) {
            minmax.push_right();
        }
        if (minmax.test() && minmax.length() > b-a+1) {
            a = minmax.left;
            b = minmax.right;
        }
        if (minmax.length() > 1)
            minmax.pop_left(); 
        else
            minmax.push_right();
    }
}

//-----------------------------------------------------------------------------
void max_admissible_subinterval(int n, int* lo, int* hi, int &a, int &b, 
                                int limit_size)
{
    if (n < limit_size)
        max_admissible_subinterval2(n, lo, hi, a, b);
    else
        max_admissible_subinterval1(n, lo, hi, a, b);
}

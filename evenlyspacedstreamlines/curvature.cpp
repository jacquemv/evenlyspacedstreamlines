#include <math.h>
#include "curvature.h"
#include "algebra.h"
#include "geometry.h"
#include "maxsubinterval.cpp"

//-----------------------------------------------------------------------------
Curvature::Curvature()
{
    points = NULL;
    constraint_lo = NULL;
    constraint_hi = NULL;
}

//-----------------------------------------------------------------------------
Curvature::Curvature(int max_nb_segments_)
{
    initialize(max_nb_segments_);
}

//-----------------------------------------------------------------------------
void Curvature::initialize(int max_nb_segments_)
{
    max_nb_segments = max_nb_segments_;
    points = new double [6*max_nb_segments];
    constraint_lo = new int [max_nb_segments];
    constraint_hi = new int [max_nb_segments];
    reset();
}

//-----------------------------------------------------------------------------
Curvature::~Curvature()
{
    delete [] points;
    delete [] constraint_lo;
    delete [] constraint_hi;
}

//-----------------------------------------------------------------------------
void Curvature::reset()
{
    nb_segments = 0;
}

//-----------------------------------------------------------------------------
void Curvature::append_segment(double* p, double* q)
{
    int k = 6*nb_segments;
    points[k  ] = p[0];
    points[k+1] = p[1];
    points[k+2] = p[2];
    points[k+3] = q[0];
    points[k+4] = q[1];
    points[k+5] = q[2];
    nb_segments++;
}

//-----------------------------------------------------------------------------
bool Curvature::calculate_constraint(double radius)
{
    bool constraint_found = 0;
    double radius2 = radius * radius;
    for (int i=0;i<nb_segments;i++) {

        double *a = points+6*i, *b = points+6*i+3; // current segment
        constraint_lo[i] = 0;
        constraint_hi[i] = nb_segments-1;

        for (int j=i+2;j<nb_segments;j++) {
            double t1, t2, *p = points+6*j, *q = points+6*j+3; // neighbor
            if (segment_distance2(p, q, a, b) > radius2)
                break;
            if (intersect_cylinder(p, q, a, b, radius, t1, t2)) {
                constraint_hi[i] = j-1;
                constraint_found = 1;
            }
        }

        for (int j=i-2;j>=0;j--) {
            double t1, t2, *p = points+6*j, *q = points+6*j+3; // neighbor
            if (segment_distance2(p, q, a, b) > radius2)
                break;
            if (intersect_cylinder(p, q, a, b, radius, t1, t2)) {
                constraint_lo[i] = j+1;
                constraint_found = 1;
            }
        }
    }
    return constraint_found;
}

//-----------------------------------------------------------------------------
void Curvature::low_curvature_interval(double radius, int &start, int &len)
{
    if (calculate_constraint(radius)) {
        int a, b;
        max_admissible_subinterval(nb_segments, constraint_lo, constraint_hi, 
                                   a, b, 100000);
        assert(is_subinterval_admissible(nb_segments, constraint_lo, 
                                         constraint_hi, a, b));
        start = a;
        len = b - a + 1;
    } else {
        start = 0;
        len = nb_segments;
    }
}

#include <math.h>
#include "quadraticforms.h"
#include "algebra.h"
#include "triangle.h"

//-----------------------------------------------------------------------------
Triangle::Triangle(double* ver_, int* idx, double* orient_, bool flip_orient,
                   bool allow_tweaking_orientation)
{
    double v1[3], v2[3], normal[3], orthogonal[3];
    int i0 = 0, i1 = 0;
    ver = ver_;

    // normal vector
    vdiff(v1, &ver[3*idx[0]], &ver[3*idx[1]]);
    vdiff(v2, &ver[3*idx[0]], &ver[3*idx[2]]);
    vcross(normal, v1, v2);
    vnormalize(normal);
    vcopy(orient, orient_);

    for (int iter=0;iter<16;iter++) {
        // orientation vector (project and normalize)
        vrmcomp(orient, normal);
        vnormalize(orient);

        // orthogonal vector
        vcross(orthogonal, normal, orient);
        if (flip_orient)
            vswap(orthogonal, orient);

        // find origin, base and apex nodes
        // i0 = argmin([0, s1, s2])
        // i1 = argmax([0, s1, s2])
        double s1, s2;
        s1 = vdot(v1, orthogonal);
        s2 = vdot(v2, orthogonal);
        i0 = 0; i1 = 0;
        if (s1 < s2) {
            if (s1 < 0) i0 = 1;
            if (s2 > 0) i1 = 2;
        } else {
            if (s2 < 0) i0 = 2;
            if (s1 > 0) i1 = 1;
        }
        int i2 = 3-i0-i1;
        // preserve orientation
        if ((i1-i0 == 2) || (i1-i0 == -1)) {
            int temp = i0;
            i0 = i1;
            i1 = temp;
        }

        // origin, base and apex
        double apex_vector[3];
        idx_origin = idx[i0];
        idx_base = idx[i1];
        idx_apex = idx[i2];
        vcopy(origin, &ver[3*idx_origin]); // origin
        vdiff(base_vector, origin, &ver[3*idx_base]); // base
        vdiff(apex_vector, origin, &ver[3*idx_apex]); // apex

        // compute s_crit and height
        double prod = vdot(orient, base_vector);
        solve2x2(prod, 1.0, vnorm2(base_vector), prod, 
                vdot(orient, apex_vector), vdot(apex_vector, base_vector), 
                s_crit, height);

        // fix sign and norm of orientation vector
        if (height != 0) vscale(orient, height);
        if (height < 0) {
            height = -height;
        }

        // precomputations
        b2 = vnorm2(base_vector);
        h2 = vnorm2(orient);
        bh = vdot(base_vector, orient);
        shrink = 1.0 / sqrt(b2 - bh*bh/h2);

        // if orientation is not parallel to an edge, everything is fine
        if (!allow_tweaking_orientation) {
            if (b2 == 0) s_crit = 0;
            break;
        }
        if (s_crit > 0 && s_crit < 1 && b2 > 0)
            break;
        
        // tweak the orientation vector and rerun
        vcopy(orient, orient_);
        double norm = sqrt(vnorm2(orient));
        if ((norm == 0) || isnan(norm)) norm = 1;
        for (int i=0;i<3;i++)
            orient[i] += norm * 1e-9*(rand()/(double)RAND_MAX-0.5);
    }

    reverse_orient = vdot(orient, orient_) < 0;
}

//-----------------------------------------------------------------------------
void Triangle::local_to_global(double s, double t, double* P)
{
    vcopyadd(P, origin, s, base_vector);
    vadd(P, t, orient);
}

//-----------------------------------------------------------------------------
void Triangle::global_to_local(double* P, double& s, double& t)
{
    double normal[3], P2[3];
    vcross(normal, base_vector, orient);
    vcopy(P2, P);
    vrmcomp(P2, normal);
    vadd(P2, -1, origin);
    solve2x2(b2, bh, bh, h2, vdot(P2, base_vector), vdot(P2, orient), s, t);
}

//-----------------------------------------------------------------------------
void Triangle::segment(double s, double* P, double* Q)
{
    vcopyadd(P, origin, s, base_vector);
    double t = (s <= s_crit) ? s/s_crit : (1-s)/(1-s_crit);
    vcopyadd(Q, P, t, orient);
}

//-----------------------------------------------------------------------------
double Triangle::segment_length(double s)
{
    if (s <= s_crit)
        return height*s/s_crit;
    else
        return height*(1-s)/(1-s_crit);
}

//-----------------------------------------------------------------------------
inline bool Triangle::contains(double s, double t)
{
    if ((s < 0) || (s > 1) || (t < 0)) return 0;
    if (s <= s_crit) return t <= s/s_crit;
    return t <= (1-s)/(1-s_crit);
}

//-----------------------------------------------------------------------------
double Triangle::vertex_coordinate(int idx)
{
    if (idx == idx_origin) return 0;
    if (idx == idx_base) return 1;
    if (idx == idx_apex) return s_crit;
    return -1;
}

//-----------------------------------------------------------------------------
void Triangle::bounding_sphere(double *center, double &radius)
{
    double *A = ver + 3*idx_origin;
    double *B = ver + 3*idx_base;
    double *C = ver + 3*idx_apex;
    double a2 = vdist2(B, C);
    double b2 = vdist2(A, C);
    double c2 = vdist2(A, B);
    // obtuse angle
    double da = b2 + c2 - a2;
    if (da <= 0) {
        radius = sqrt(a2)/2;
        vlincomb(center, 0.5, B, 0.5, C);
        return;
    }
    double db = a2 + c2 - b2;
    if (db <= 0) {
        radius = sqrt(b2)/2;
        vlincomb(center, 0.5, A, 0.5, C);
        return;
    }
    double dc = a2 + b2 - c2;
    if (dc <= 0) {
        radius = sqrt(c2)/2;
        vlincomb(center, 0.5, A, 0.5, B);
        return;
    }
    // circumscribed circle
    vlincomb(center, a2*da, A, b2*db, B);
    vadd(center, c2*dc, C);
    vscale(center, 1.0/(a2*da+b2*db+c2*dc));
    radius = vdist(center, A);
}

//-----------------------------------------------------------------------------
int Triangle::common_edge(int* idx)
{
    bool O, A, B;
    O = (idx[0] == idx_origin) || (idx[1] == idx_origin) 
                               || (idx[2] == idx_origin);
    B = (idx[0] == idx_base) || (idx[1] == idx_base) || (idx[2] == idx_base);
    A = (idx[0] == idx_apex) || (idx[1] == idx_apex) || (idx[2] == idx_apex);
    if (O && B) return 0;
    if (O && A) return 1;
    if (A && B) return 2;
    return -1;
}

//-----------------------------------------------------------------------------
int Triangle::intersection_segment_sphere(double* P, double R,
                                 double s1, double t1, double s2, double t2,
                                 double& si1, double& ti1, 
                                 double& si2, double& ti2)
{
    double x0[3];
    vdiff(x0, P, origin);
    EllipticQuadraticForm qf(b2, h2, bh,
                             vdot(base_vector, x0), vdot(orient, x0), 
                             vnorm2(x0)-R*R);
    return qf.intersection_segment(s1, t1, s2, t2, si1, ti1, si2, ti2);
}

//-----------------------------------------------------------------------------
int Triangle::coordinate_mask_quadratics(EllipticQuadraticForm qf,
                                         double& s_min, double& s_max)
{
    s_min = 9e9; s_max = -9e9;
    #define PUSH(s) \
        { if (s < s_min) s_min = s; if (s > s_max) s_max = s; }

    #define UPDATE_INTERVAL(n, s1, s2) \
        { if (n > 0) PUSH(s1); if (n > 1) PUSH(s2); }
    
    // check the vertices of the triangle
    if (qf.evaluate(1, 0) <= 0) PUSH(1);
    if (qf.evaluate(0, 0) <= 0) PUSH(0);

    // intersections with the edges of the triangle
    int n;
    double s1, t1, s2, t2;
    n = qf.intersection_segment(0, 0, 1, 0, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, s2);
    n = qf.intersection_segment(0, 0, s_crit, 1, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, s2);
    n = qf.intersection_segment(1, 0, s_crit, 1, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, s2);

    // the extrema may be in the interior of the triangle
    if (qf.extreme_points_xaxis(s1, t1, s2, t2) == 2) {
        if (contains(s1, t1)) PUSH(s1);
        if (contains(s2, t2)) PUSH(s2);
    }

    #undef PUSH
    #undef UPDATE_INTERVAL

    return 0;
}

//-----------------------------------------------------------------------------
int Triangle::coordinate_mask_quadratics_constrained(EllipticQuadraticForm qf,
                                double* line, double& s_min, double& s_max)
{
    s_min = 9e9; s_max = -9e9;
    #define PUSH(s, t) \
        { double w = line[0]*s+line[1]*t+line[2]; \
          if ((w > 0) && (w < 1)) { \
            if (s < s_min) s_min = s; \
            if (s > s_max) s_max = s; \
          } \
        }

    #define UPDATE_INTERVAL(n, s1, t1, s2, t2) \
        { if (n > 0) PUSH(s1, t1); if (n > 1) PUSH(s2, t2); }
    
    // check the vertices of the triangle
    if (qf.evaluate(1, 0) <= 0) PUSH(1, 0);
    if (qf.evaluate(0, 0) <= 0) PUSH(0, 0);

    // intersections with the edges of the triangle
    int n;
    double s1, t1, s2, t2;
    n = qf.intersection_segment(0, 0, 1, 0, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, t1, s2, t2);
    n = qf.intersection_segment(0, 0, s_crit, 1, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, t1, s2, t2);
    n = qf.intersection_segment(1, 0, s_crit, 1, s1, t1, s2, t2);
    UPDATE_INTERVAL(n, s1, t1, s2, t2);

    // the extrema may be in the interior of the triangle
    if (qf.extreme_points_xaxis(s1, t1, s2, t2) == 2) {
        if (contains(s1, t1)) PUSH(s1, t1);
        if (contains(s2, t2)) PUSH(s2, t2);
    }

    #undef PUSH
    #undef UPDATE_INTERVAL

    return 0;
}

//-----------------------------------------------------------------------------
void Triangle::linear_constraint(double* P, double* Q, double* line)
{
    double d[3];
    vdiff(d, P, Q);
    vscale(d, 1.0/vnorm2(d));
    double x0[3];
    vdiff(x0, P, origin);
    line[0] = vdot(base_vector, d);
    line[1] = vdot(orient, d);
    line[2] = vdot(x0, d);
}

//-----------------------------------------------------------------------------
int Triangle::coordinate_mask_sphere(double* P, double R, 
                                     double& s_min, double& s_max)
{
    double x0[3];
    vdiff(x0, P, origin);
    EllipticQuadraticForm qf(b2, h2, bh,
                             vdot(base_vector, x0), vdot(orient, x0), 
                             vnorm2(x0)-R*R);
    return coordinate_mask_quadratics(qf, s_min, s_max);
}

//-----------------------------------------------------------------------------
int Triangle::coordinate_mask_infinite_cylinder(double* P, double* Q, double R, 
                                       double& s_min, double& s_max)
{
    double d[3];
    vdiff(d, P, Q);
    vnormalize(d);
    double x0[3];
    vdiff(x0, P, origin);
    vadd(x0, -vdot(x0, d), d);
    double bd[3], hd[3];
    vcopyadd(bd, base_vector, -vdot(base_vector, d), d);
    vcopyadd(hd, orient, -vdot(orient, d), d);
    double b2d = vdot(bd, bd);
    double h2d = vdot(hd, hd);
    double bhd = vdot(bd, hd);
    EllipticQuadraticForm qf(b2d, h2d, bhd,
                             vdot(bd, x0), vdot(hd, x0), 
                             vnorm2(x0)-R*R);
    return coordinate_mask_quadratics(qf, s_min, s_max);
}

//-----------------------------------------------------------------------------
int Triangle::coordinate_mask_cylinder(double* P, double* Q, double R, 
                                       double& s_min, double& s_max)
{
    double d[3];
    vdiff(d, P, Q);
    vnormalize(d);
    double x0[3];
    vdiff(x0, P, origin);
    vadd(x0, -vdot(x0, d), d);
    double bd[3], hd[3];
    vcopyadd(bd, base_vector, -vdot(base_vector, d), d);
    vcopyadd(hd, orient, -vdot(orient, d), d);
    double b2d = vdot(bd, bd);
    double h2d = vdot(hd, hd);
    double bhd = vdot(bd, hd);
    EllipticQuadraticForm qf(b2d, h2d, bhd,
                             vdot(bd, x0), vdot(hd, x0), 
                             vnorm2(x0)-R*R);
    double line[3];
    linear_constraint(P, Q, line);
    coordinate_mask_quadratics_constrained(qf, line, s_min, s_max);

    double s1, s2;
    coordinate_mask_sphere(P, R, s1, s2);

    if (s1 < s_min) s_min = s1; 
    if (s2 > s_max) s_max = s2;

    coordinate_mask_sphere(Q, R, s1, s2);
    if (s1 < s_min) s_min = s1; 
    if (s2 > s_max) s_max = s2;
    return 0;
}
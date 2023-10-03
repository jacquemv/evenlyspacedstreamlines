#include "engine.cpp"

const double EPS = 1e-10;

bool allclose(double* P, double* Q)
{
    return vdist(P, Q) < EPS;
}

bool is_circshift(int a1, int a2, int a3, int b1, int b2, int b3)
{
    return ((a1 == b1) && (a2 == b2) && (a3 == b3))
        || ((a1 == b2) && (a2 == b3) && (a3 == b1))
        || ((a1 == b3) && (a2 == b1) && (a3 == b2));
}

void test_segments()
{
    double ver[] = {
        9e9, 9e9, 9e9,
        1, 0.5, 2,
        5, 1.5, 4,
        2, 4, 3
    };
    int tri[] = {
        1, 2, 3,
        2, 3, 1,
        3, 1, 2,
        3, 2, 1
    };
    double orient[] = {
        3.7 , -0.05,  1.7
    };
    for (int k=0;k<2;k++) {
        for (int i=0;i<4;i++) {
            Triangle T(ver, tri + 3*i, orient, 0, 0);
            double P[3], Q[3];

            assert( is_circshift(tri[3*i], tri[3*i+1], tri[3*i+2], T.idx_origin, T.idx_base, T.idx_apex));

            T.local_to_global(0, 0, P);
            assert( allclose(T.origin, T.ver + 3*T.idx_origin) );
            T.local_to_global(1, 0, P);
            assert( allclose(P, T.ver + 3*T.idx_base) );
            T.local_to_global(T.s_crit, 1, P);
            assert( allclose(P, T.ver + 3*T.idx_apex) );

            T.segment(T.s_crit, P, Q);
            assert( allclose(Q, T.ver + 3*T.idx_apex) );
            assert( fabs(vdist(P, Q) - T.height) < EPS );

            T.segment(0.3, P, Q);
            assert( fabs(vdist(P, Q) - T.segment_length(0.3)) < EPS );
        }
        // mirror
        for (int j=0;j<12;j++) ver[j] = -ver[j];
    }
}

void test_contains()
{
    double ver[] = {
        1, 0.5, 2,
        5, 1.5, 4,
        2, 4, 3
    };
    int tri[] = {0, 1, 2};
    double orient[] = {3.7 , -0.05,  1.7};
    Triangle T(ver, tri, orient, 0, 0);
    assert( T.contains(1, 0) );
    assert( T.contains(0, 0) );
    assert( ! T.contains(-0.1, 0) );
    assert( ! T.contains(2, 0) );
    assert( ! T.contains(0.5, -1) );
    assert( ! T.contains(0.2, 1.01) );
    assert( T.contains( T.s_crit, 0.999) );
    assert( T.contains( T.s_crit/2, 0.2) );
    assert( ! T.contains( T.s_crit/2, 0.55) );
    assert( T.contains( T.s_crit + 0.5*(1-T.s_crit), 0.3) );
    assert( ! T.contains( T.s_crit + 0.5*(1-T.s_crit), 0.7) );
}

void test_projection()
{
    double ver[] = {
        1, 0.5, 2,
        5, 1.5, 4,
        2, 4, 3
    };
    int tri[] = {0, 1, 2};
    double orient[] = {3.7 , -0.05,  1.7};
    Triangle T(ver, tri, orient, 0, 0);
    double P[3], s, t;
    T.local_to_global(0.2, 0.5, P);
    T.global_to_local(P, s, t);
    assert( fabs(s-0.2)+fabs(t-0.5) < EPS );
}

void test_intersect_sphere()
{
    double ver[] = {
        1, 0.5, 2,
        5, 1.5, 4,
        2, 4, 3
    };
    int tri[] = {0, 1, 2};
    double orient[] = {3.7 , -0.05,  1.7};
    Triangle T(ver, tri, orient, 0, 0);
    double center[] = {5, 2, 5};
    double s1, t1, s2, t2;
    double P1[3]; T.local_to_global(0.1, -0.2, P1);
    double P2[3]; T.local_to_global(1.2, 2.5, P2);
    //printf("%.3f %.3f %.3f\n", P1[0], P1[1], P1[2]);
    //printf("%.3f %.3f %.3f\n", P2[0], P2[1], P2[2]);
    double R = 1.3;
    int n = T.intersection_segment_sphere(center, R, 0.1, -0.2, 1.2, 2.5, s1, t1, s2, t2);
    assert( n == 2 );
    if (n >= 1) {
        double P[3];
        T.local_to_global(s1, t1, P);
        assert( fabs(vdist(center, P)-R) < EPS);
    }
    if (n >= 2) {
        double P[3];
        T.local_to_global(s2, t2, P);
        assert( fabs(vdist(center, P)-R) < EPS);
    }
}

void test_mask_sphere()
{
    double ver[] = {
        2, 1, 0,
        5, 1, 0,
        3, 3, 0
    };
    int tri[] = {0, 1, 2};
    double orient[] = {0, 1, 0};
    Triangle T(ver, tri, orient, 0, 0);
    double P[3] = {3, 2, 0};
    double s1, s2;
    #define TEST_SPHERE(R, xmin, xmax) \
        T.coordinate_mask_sphere(P, R, s1, s2); \
        if (s2 < s1) \
            assert( xmax < xmin ); \
        else \
            assert( fabs(2+3*s1 - xmin) + fabs(2+3*s2 - xmax) < EPS );
    TEST_SPHERE(0.2, 2.8, 3.2);
    TEST_SPHERE(0.5, 2.5, 3.5);
    TEST_SPHERE(1, 2.2, 4);
    TEST_SPHERE(1.5, 2, 4.435414346693484);
    TEST_SPHERE(3, 2, 5);
    P[1] = 0;
    TEST_SPHERE(1, 1, -1);
    TEST_SPHERE(2, 2, 4.732050807568877);
    TEST_SPHERE(2.5, 2, 5);
    P[0] = 4;
    TEST_SPHERE(1.25, 3.25, 4.75);

    //T.coordinate_mask_sphere(P, 1.25, s1, s2);
    //printf("%.15f %.15f\n", 2+3*s1, 2+3*s2);
}

void test_mask_cylinder()
{
    double ver[] = {
        2, 1, 0,
        5, 1, 0,
        3, 3, 0
    };
    int tri[] = {0, 1, 2};
    double orient[] = {0, 1, 0};
    Triangle T(ver, tri, orient, 0, 0);
    //double P[3] = {2, 0, 0}, Q[3] = {5, 3, 0};
    double s1, s2;
    #define TEST_CYLINDER(R, xmin, xmax) \
        T.coordinate_mask_cylinder(P, Q, R, s1, s2); \
        if (s2 < s1) \
            assert( xmax < xmin ); \
        else \
            assert( fabs(2+3*s1 - xmin) + fabs(2+3*s2 - xmax) < EPS );
    
    double R = 0.25;
    double P[3] = {2, 2, 0}, Q[3] = {4.5, 2, 0};
    TEST_CYLINDER(0.25, 2.375, 4.25);
    P[1] = Q[1] = 2.5;
    TEST_CYLINDER(0.75, 2.375, 4.25);
    P[1] = Q[1] = 0.75;
    TEST_CYLINDER(1, 2, 5);
    P[0] = 2; P[1] = 0; Q[0] = 5; Q[1] = 3;
    TEST_CYLINDER(0.5, 2.292893218813453, 4.353553390593274);

    P[0] = 2.5; P[1] = 1.5; Q[0] = 3.5; Q[1] = 2;
    TEST_CYLINDER(0.2, 2.3, 3.7);
    TEST_CYLINDER(0.25, 2.25, 3.75);
    TEST_CYLINDER(0.6, 2.051002008040225, 4.092782730020049);
    
    P[0] = 2; P[1] = 0; Q[0] = 3.5; Q[1] = 1.5;
    TEST_CYLINDER(0.5, 2.292893218813453, 4);

    //R = 0.5;
    //T.coordinate_mask_cylinder(P, Q, R, s1, s2);
    //printf("%.15f %.15f\n", 2+3*s1, 2+3*s2);
}

void test_linear_constraint()
{
    double ver[] = {
        2, 1, 0,
        5, 1, 0,
        3, 3, 0
    };
    int tri[] = {0, 1, 2};
    double orient[] = {0, 1, 0};
    Triangle T(ver, tri, orient, 0, 0);
    double P[3] = {3, 0.5, 0}, Q[3] = {4, 1.5, 0};
    double line[3];
    T.linear_constraint(P, Q, line);
    //printf("%.5f, %.5f, %.5f\n", line[0], line[1], line[2]);
}

void test_shrink()
{
    double ver[] = {
        1, 0.5, 2,
        5, 1.5, 4,
        2, 4, 3
    };
    int tri[] = {0, 1, 2};
    double orient[] = {3.7 , -0.05,  1.7};
    Triangle T(ver, tri, orient, 0, 0);
    double P1[3], Q1[3], P2[3], Q2[3];
    T.segment(0.6, P1, Q1);
    T.segment(0.61, P2, Q2);
    double d = sqrt(segment_distance2(P1, Q1, P2, Q2));
    //printf("%.5f, %.5f\n", d, T.shrink);
    assert( fabs(0.01/d - T.shrink) < EPS );
}

void test_enclosing_circle()
{
    const double EPS = 1e-6;
    double ver[] = {
        2, 1, 1,
        5, 1, 2,
        3, 3, -1
    };
    int tri[] = {0, 1, 2};
    double orient[] = {0, 1, 0};
    Triangle T(ver, tri, orient, 0, 0);
    double c[3], r;
    T.bounding_sphere(c, r);
    assert( fabs(c[0] - 3.81460674) < EPS );
    assert( fabs(c[1] - 1.8988764) < EPS );
    assert( fabs(c[2] - 0.55617978) < EPS );
    assert( fabs(r - 2.07310222) < EPS );
}

void generate_random_enclosing_circles()
{
    #define RND (double)(rand())/RAND_MAX
    for (int i=0;i<100;i++) {
        double ver[] = {
            RND, RND, 0,
            RND, RND, 0,
            RND, RND, 0
        };
        int tri[] = {0, 1, 2};
        double orient[] = {0, 1, 0};
        Triangle T(ver, tri, orient, 0, 0);
        double c[3], r;
        T.bounding_sphere(c, r);
        printf("%.3f %.3f   %.3f %.3f  %.3f %.3f   ", ver[0], ver[1], ver[3], ver[4], ver[6], ver[7]);
        printf("%.3f %.3f %.3f  ", c[0], c[1], c[2]);
        printf("%.3f\n", r);
    }
}

int main()
{
    test_segments();
    test_contains();
    test_projection();
    test_intersect_sphere();
    test_mask_sphere();
    test_mask_cylinder();
    test_linear_constraint();
    test_shrink();
    test_enclosing_circle();
    //generate_random_enclosing_circles();
}
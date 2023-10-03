#include "engine.cpp"

int main()
{
    const double EPS = 1e-10;
    double param[][5] = {
        // xc, yc, a, b, th
        {-2.3, 1.2, 2.5, 0.7, 1.1},
        {-2.3, 1.2, 2.5, 0.7, -0.3},
        {1, 0, 250, 0.7, 0},
        {2.3, -1.2, 200, 0.1, 2*atan(1)},
        {0, 0, 1, 1, 0},
    };

    // ellipses
    for (int i=0;i<5;i++) {
        double xc = param[i][0], yc = param[i][1];
        double a = param[i][2], b = param[i][3], th = param[i][4];
        EllipticQuadraticForm qf;
        qf.ellipse(xc, yc, a, b, th);
        double x, y, x1, y1, x2, y2;
        assert( qf.delta > 0 );
        // test parameters
        assert( qf.center(x, y) == 1 );
        assert( (fabs(x-xc)<EPS) && (fabs(y-yc)<EPS) );
        assert( qf.semiaxis(x, y) == 1 );
        assert( (fabs(x-a)<EPS) && (fabs(y-b)<EPS) );
        assert( qf.halfwidth(x) == 1 );
        assert( fabs(x-b)<EPS );
        assert( qf.angle(x) == 1 );
        if (a != b)
            assert( fabs(x-th)<EPS );
        // test minimum
        qf.center(x, y);
        assert( fabs(qf.evaluate(x, y) - qf.minimum()) < EPS );
        // extreme points x axis
        assert( qf.extreme_points_xaxis(x1, y1, x2, y2) == 2);
        assert( x1 <= x2 );
        assert( fabs(qf.evaluate(x1, y1)) < EPS );
        assert( fabs(qf.evaluate(x2, y2)) < EPS );
        // extreme points y axis
        assert( qf.extreme_points_yaxis(x1, y1, x2, y2) == 2);
        assert( y1 <= y2 );
        assert( fabs(qf.evaluate(x1, y1)) < EPS );
        assert( fabs(qf.evaluate(x2, y2)) < EPS );
    }

    // empty set
    for (int i=0;i<5;i++) {
        double xc = param[i][0], yc = param[i][1];
        double a = param[i][2], b = param[i][3], th = param[i][4];
        EllipticQuadraticForm qf_base;
        double x, y, x1, y1, x2, y2;
        qf_base.ellipse(xc, yc, a, b, th);
        assert( qf_base.determinant < 0 );
        EllipticQuadraticForm qf(qf_base.a, qf_base.b, qf_base.c, 
                                 qf_base.d, qf_base.e, 
                                 qf_base.f - qf_base.minimum() + 1e-6);
        assert( qf.determinant > 0 );
        assert( qf.center(x, y) == 1 );
        assert( (fabs(x-xc)<EPS) && (fabs(y-yc)<EPS) );
        assert( fabs(qf.evaluate(x, y) - qf.minimum()) < EPS );
        assert( (fabs(x-xc)<EPS) && (fabs(y-yc)<EPS) );
        assert( qf.angle(x) == 1 );
        if (a != b)
            assert( fabs(x-th)<EPS );
        assert( qf.semiaxis(x, y) == 0 );
        assert( qf.extreme_points_xaxis(x1, y1, x2, y2) == 0);
        assert( qf.extreme_points_yaxis(x1, y1, x2, y2) == 0);
    }

    // infinite strip
    double paramb[][4] = {
        // xc, yc, w, th
        {-2.3, 1.2, 2.5, 0.7},
        {-0.6, 0.3, 1, 0},
        {0.7, 0.2, 3, 2*atan(1)},
    };
    for (int i=0;i<3;i++) {
        double xc = paramb[i][0], yc = paramb[i][1];
        double w = paramb[i][2], th = paramb[i][3];
        EllipticQuadraticForm qf;
        qf.infinite_strip(xc, yc, w, th);
        assert( qf.delta == 0 );
        assert( qf.determinant == 0 );
        double x, y, x1, y1, x2, y2;
        assert( qf.center(x, y) == 1 );
        assert( fabs(qf.evaluate(x, y) - qf.minimum()) < EPS );
        assert( qf.angle(x) == 1 );
        assert( fabs(x-th)<EPS );
        assert( qf.semiaxis(x, y) == 0 );
        assert( qf.halfwidth(x) == 1 );
        assert( fabs(x-w)<EPS );
        assert( qf.extreme_points_xaxis(x1, y1, x2, y2) == 0);
        assert( qf.extreme_points_yaxis(x1, y1, x2, y2) == 0);
    }

    // empty set
    for (int i=0;i<3;i++) {
        double xc = paramb[i][0], yc = paramb[i][1];
        double w = paramb[i][2], th = paramb[i][3];
        EllipticQuadraticForm qf_base;
        qf_base.infinite_strip(xc, yc, w, th);
        EllipticQuadraticForm qf(qf_base.a, qf_base.b, qf_base.c, 
                                 qf_base.d, qf_base.e, 
                                 qf_base.f - qf_base.minimum() + 1e-6);
        assert( qf.delta == 0 );
        assert( fabs(qf.determinant) < EPS );
        double x, y, x1, y1, x2, y2;
        assert( qf.center(x, y) == 1 );
        assert( fabs(qf.evaluate(x, y) - qf.minimum()) < EPS );
        assert( qf.angle(x) == 1 );
        assert( fabs(x-th)<EPS );
        assert( qf.semiaxis(x, y) == 0 );
        assert( qf.halfwidth(x) == 0 );
        assert( qf.extreme_points_xaxis(x1, y1, x2, y2) == 0);
        assert( qf.extreme_points_yaxis(x1, y1, x2, y2) == 0);
    }

    // intersections
    EllipticQuadraticForm qf;
    qf.ellipse(-0.2, 0.7, 2, 0.4, -2);
    double segments[][5] = {
        {-0.5, 0, -0.5, 0.75,  0},
        {-0.5, -1, -0.5, 0.75,  1},
        {-0.5, 0.75, -0.5, -1,  1},
        {-0.5, -0.3, -0.5, 2,  1},
        {-0.5, -4, -0.5, 5,  2},
        {0.25, 1, 3, 1,  1},
        {3, 1, 0.25, 1,  1},
        {0.25, -0.5, -1, 0.1,  2},
        {0.25, -0.5, -1, 0.1,  2},
        {0.5, 0.6, 0.7, 0.9,  0},
    };
    for (int i=0;i<8;i++) {
        double x, y, x1, y1, x2, y2;
        int n;
        n = qf.intersection_segment(segments[i][0], segments[i][1], 
                                    segments[i][2], segments[i][3],
                                    x1, y1, x2, y2);
        assert( n == segments[i][4] );
        if (n > 0) assert( fabs(qf.evaluate(x1, y1)) < EPS );
        if (n > 1) assert( fabs(qf.evaluate(x2, y2)) < EPS );
        // swap points
        n = qf.intersection_segment(segments[i][0], segments[i][1], 
                                    segments[i][2], segments[i][3],
                                    x2, y2, x1, y1);
        assert( n == segments[i][4] );
        if (n > 0) assert( fabs(qf.evaluate(x1, y1)) < EPS );
        if (n > 1) assert( fabs(qf.evaluate(x2, y2)) < EPS );
    }
    // infinite strip
    double params[][13] = {
        // xc, yc, w, th,    x1, y1, x2, y2,    n
        {0, 0, 1, atan(1),   2.41, 1, 2.4145, 1,  1},
        {0, 0, 1, atan(1),   2, 1, 2.41, 1,  0},
        {0, 0, 1, atan(1),   -5, 0, 5, 2,  2},
        {0, 0, 1, atan(1),   3, 4.41, 3.000001, 4.4145,  1},
        {0, 0, 1, atan(1),   3, 4.41, 3, 4.4145,  1},
        {0, 0, 1, atan(1),   3, 4.41, 3, 4.4145,  1},
        {0, 0, 1, atan(1),   3, -4.41, 3, 3,  1},
        {0, 0, 1, atan(1),   3, 4.46, 3, -4.4145,  2},
        {0, 0, 1, 2*atan(1),   0, -100, 0, 100,  0},
        {0.5, 0.5, 2, 0,   1, 1, 3, 3,  1},
        {0.5, 0.5, 2, 0,   3, -1, 1, -2,  1},
        {0.5, 0.5, 2, 0,   1, -3, 1, 4,  2},
        {0.5, 0.5, 2, 0,   2, -3, 1, 4,  2},
    };
    for (int i=0;i<10;i++) {
        double x, y, x1, y1, x2, y2;
        int n;
        double xc = params[i][0], yc = params[i][1];
        double w = params[i][2], th = params[i][3];
        EllipticQuadraticForm qf;
        qf.infinite_strip(xc, yc, w, th);
        if (i==8) qf.c = 0; // fix inaccuracy of sin/cos at pi/2
        n = qf.intersection_segment(params[i][4], params[i][5], 
            params[i][6], params[i][7], x1, y1, x2, y2);
        assert( n == params[i][8] );
        if (n > 0) assert( fabs(qf.evaluate(x1, y1)) < EPS );
        if (n > 1) assert( fabs(qf.evaluate(x2, y2)) < EPS );
        // swap points
        n = qf.intersection_segment(params[i][4], params[i][5], 
            params[i][6], params[i][7], x2, y2, x1, y1);
        assert( n == params[i][8] );
        if (n > 0) assert( fabs(qf.evaluate(x1, y1)) < EPS );
        if (n > 1) assert( fabs(qf.evaluate(x2, y2)) < EPS );
    }
}
#include "engine.cpp"


void curv_calc(int nb_points, double* points, double radius, 
               int &start, int &len, bool print_intervals, bool print_points)
{
    Curvature curv(nb_points);
    curv.reset();
    for (int i=0;i<nb_points-1;i++)
        curv.append_segment(points+3*i, points+3*i+3);
    curv.low_curvature_interval(radius, start, len);

    if (print_points) {
        printf("%d %.2f\n", nb_points, radius);
        for (int i=0;i<nb_points;i++)
            printf("%.4f %.4f %.4f\n", points[3*i], points[3*i+1], points[3*i+2]);
        printf("%d %d\n", start, len);
    }

    if (print_intervals) {
        for (int i=0;i<curv.nb_segments;i++) {
            printf("%d: %d %d\n", i, curv.constraint_lo[i], curv.constraint_hi[i]);
        }
        printf("start = %d, len = %d, interval = [%d, %d]\n", start, len, start, start+len-1);
    }
}

void test_basic()
{
    double points[] = {
        -1, 1, 0,
        0, 0, 0,
        2, 0, 0,
        3, 1, 0,
        2.5, 2, 0,
        1, 1.5, 0
    };
    int start, len;
    curv_calc(6, points, 2, start, len, 0, 0);
    assert(start == 0 && len == 4);
}

void test_ellipse()
{
    int n = 20;
    double *points = new double [3*n];
    for (int i=0;i<n;i++) {
        points[3*i]   = cos(2.0*M_PI/n*i);
        points[3*i+1] = 2*sin(2.0*M_PI/n*i);
        points[3*i+2] = 0;
    }

    int start, len;
    curv_calc(n, points, 1.7, start, len, 0, 1);
    delete [] points;
}

void test_U_turn()
{
    int m = 10, n = 3*m;
    double *points = new double [3*n];
    for (int i=0;i<m;i++) {
        points[3*(i)]   = 2*(i/((double)m)-1);
        points[3*(i)+1] = 1;
        points[3*(i)+2] = 0;
    }
    for (int i=0;i<m;i++) {
        points[3*(i+m)]   = sin(M_PI/m*i);
        points[3*(i+m)+1] = cos(M_PI/m*i);
        points[3*(i+m)+2] = 0;
    }
    for (int i=0;i<m;i++) {
        points[3*(i+2*m)]   = 2*(-i/((double)m));
        points[3*(i+2*m)+1] = -1;
        points[3*(i+2*m)+2] = 0;
    }

    int start, len;
    curv_calc(n, points, 3.0, start, len, 0, 1);
    delete [] points;
}

int main(int argc, char** argv)
{
    test_basic();
    test_U_turn();
}
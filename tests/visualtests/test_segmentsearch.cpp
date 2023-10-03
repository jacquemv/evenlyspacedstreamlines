#include "algebra.cpp"
#include "geometry.cpp"
#include "segmentsearch.cpp"
#include <math.h>
#include <stdio.h>
#include <vector>

int main(int argc, char** argv)
{
    double b = atof(argv[1]); // 0.95
    double step = atof(argv[2]); // 0.05;
    double radius = atof(argv[3]); // 0.05;
    int N = atoi(argv[4]); // 1000

    double P[3], Q[3];
    P[2] = Q[2] = 0;
    P[0] = 1;
    P[1] = 0;
    printf("%.4f %.4f\n", P[0], P[1]);

    int nb_neigh[] = {2, 2, 2};
    int neigh[] = {1, 2, 999, 999, 
                   0, 2, 999, 999, 
                   0, 1, 999, 999};
    int* neigh_list[] = {neigh, neigh+4, neigh+8};
    SegmentSearch collection;
    collection.initialize(3, nb_neigh, neigh_list);
    
    for (int i=0;i<N;i++) {
        double theta = i*step;
        double r = pow(b, theta);
        Q[0] = r * cos(theta);
        Q[1] = r * sin(theta);
        if (collection.is_too_close(P, Q, i % 3, radius))
            break;
        collection.insert(P, Q, i % 3);
        printf("%.4f %.4f\n", Q[0], Q[1]);
        P[0] = Q[0];
        P[1] = Q[1];
    }

    collection.reverse_order();
    printf("9e9 9e9\n");

    P[0] = 1;
    P[1] = 0;
    for (int i=0;i<N;i++) {
        double theta = i*step;
        double r = pow(b, theta);
        Q[0] = 2 - r * cos(theta);
        Q[1] = - r * sin(theta);
        if (collection.is_too_close(P, Q, i % 3, radius))
            break;
        collection.insert(P, Q, i % 3);
        printf("%.4f %.4f\n", Q[0], Q[1]);
        P[0] = Q[0];
        P[1] = Q[1];
    }
}
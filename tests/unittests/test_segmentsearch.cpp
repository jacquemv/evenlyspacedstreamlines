#include "engine.cpp"
#include <vector>
#include <cassert>

void test_query()
{
    const int max_segments = 10;
    int nb_bins = 5;
    int bin_maxsize = 4;
    int nb_neigh[] = {0, 2, 2, 3, 0};
    int neigh[] = {999, 999, 999,   0, 2, 999,   1, 4, 999,   1, 2, 4,   999, 999, 999};
    int* neigh_list[] = {neigh, neigh+3, neigh+6, neigh+9, neigh+12};
    SegmentSearch collection;
    collection.initialize(nb_bins, nb_neigh, neigh_list);

    double points[6*max_segments];
    for (int i=0;i<6*max_segments;i++) points[i] = i;

    #define INSERT(bin_idx, pt_idx) collection.insert(points+6*pt_idx, points+6*pt_idx+3, bin_idx);
    INSERT(0, 1);
    INSERT(2, 5); INSERT(2, 6);
    INSERT(1, 3);
    INSERT(4, 0); INSERT(4, 4); INSERT(4, 8);
    INSERT(3, 9); INSERT(3, 2); INSERT(3, 7);

    std::vector<SegmentNode> seglist;
    for (int bin_idx=0;bin_idx<nb_bins;bin_idx++) {
        printf("bin %d:\n", bin_idx);
        collection.query_neighborhood(bin_idx, seglist);
        for (unsigned int i=0;i<seglist.size();i++) {
            double *p = seglist[i].P, *q = seglist[i].Q;
            printf("%.0f %.0f %.0f   %.0f %.0f %.0f\n", p[0], p[1], p[2], q[0], q[1], q[2]);
        }
    }
}

void test_too_close_within_a_bin()
{
    int nb_bins = 5;
    int bin_maxsize = 4;
    int nb_neigh[] = {0, 2, 2, 3, 0};
    int neigh[] = {999, 999, 999,   0, 2, 999,   1, 4, 999,   1, 2, 4,   999, 999, 999};
    int* neigh_list[] = {neigh, neigh+3, neigh+6, neigh+9, neigh+12};
    SegmentSearch collection;
    collection.initialize(nb_bins, nb_neigh, neigh_list);

    int idx = 3;
    double X[] = {0, 1, 0,      1, 1, 0,
                  0, 1.5, 0,    1, 1.5, 0,
                  0, 1.8, 0,    1, 1.8, 0,
                  0, 1.65, 0,   1, 1.65, 0,
                  0, 1.76, 0,   1, 1.76, 0,
                  0, 3, 0,      1, 3, 0,
                  0, 1.82, 0,   1, 1.82, 0};
    double r = 0.1;
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X, X+3, idx) == 1);
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X+6, X+9, idx) == 2);
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X+12, X+15, idx) == 3);
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X+18, X+21, idx) == 4);
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X+24, X+27, idx) == 5);
    assert( collection.too_close_within_a_bin(idx, r) == 1 );
    assert( collection.insert(X+30, X+33, idx) == 6);
    assert( collection.too_close_within_a_bin(idx, r) == 0 );
    assert( collection.insert(X+36, X+39, idx) == 7);
    assert( collection.too_close_within_a_bin(idx, r) == 1 );

}

int main()
{
    test_too_close_within_a_bin();
    test_query();
}
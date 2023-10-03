#ifndef SEGMENTSEARCH_H_
#define SEGMENTSEARCH_H_

#include "geometry.h"
#include <vector>


//-----------------------------------------------------------------------------
typedef struct SegmentNode {
    int id;
    double P[3], Q[3];
} SegmentNode;

//-----------------------------------------------------------------------------
class SegmentSearch {
public:
    SegmentSearch();
    ~SegmentSearch();
    void initialize(int nb_bins_, int* bin_nb_neigh_, int** bin_neigh_);
    void preallocate(int bin_size);
    void clear();
    void reverse_order();

    int insert(double* P, double* Q, int bin_idx);
    void query(int bin_idx, std::vector<SegmentNode> &seg_list);
    void query_neighborhood(int bin_idx, std::vector<SegmentNode> &seg_list);

    bool is_too_close(double* P, double* Q, int bin_idx, double radius);
    void find_too_close(double* P, double* Q, int bin_idx, 
                        double radius, std::vector<SegmentNode> &seg_list);
    bool too_close_within_a_bin(int idx, double radius);

//private:
    // create 'nb_bins bins'
    int nb_bins;
    std::vector<SegmentNode> *bin;

    // number of neighbors for each bin
    int* bin_nb_neigh; 
    // nb_bins arrays defining the neighborhood of each bin
    int** bin_neigh;

    int next_id; // segment id
    std::vector<int> selection;
};

#endif
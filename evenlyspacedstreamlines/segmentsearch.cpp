#include <stdlib.h>
#include "geometry.h"
#include "segmentsearch.h"


//-----------------------------------------------------------------------------
SegmentSearch::SegmentSearch()
{
    nb_bins = 0;
    bin = (std::vector<SegmentNode> *) NULL;
}

//-----------------------------------------------------------------------------
SegmentSearch::~SegmentSearch()
{
    delete [] bin;
}

//-----------------------------------------------------------------------------
void SegmentSearch::initialize(int nb_bins_, int* bin_nb_neigh_, 
                               int** bin_neigh_)
{
    nb_bins = nb_bins_;
    bin_nb_neigh = bin_nb_neigh_;
    bin_neigh = bin_neigh_;
    bin = new std::vector<SegmentNode> [nb_bins];
    next_id = 0;
}

//-----------------------------------------------------------------------------
void SegmentSearch::preallocate(int bin_size)
{
    for (int i=0;i<nb_bins;i++) {
        bin[i].reserve(bin_size);
    }
}

//-----------------------------------------------------------------------------
void SegmentSearch::clear()
{
    for (int i=0;i<nb_bins;i++) {
        bin[i].clear();
    }
}

//-----------------------------------------------------------------------------
void SegmentSearch::reverse_order()
{
    for (int i=0;i<nb_bins;i++)
        for (size_t k=0;k<bin[i].size();k++) {
            SegmentNode* s = &bin[i][k];
            s->id = next_id-1 - s->id;
        }
}

//-----------------------------------------------------------------------------
int SegmentSearch::insert(double* P, double* Q, int bin_idx)
{
    SegmentNode item = {next_id, {P[0], P[1], P[2]}, {Q[0], Q[1], Q[2]}};
    bin[bin_idx].push_back(item);
    next_id++;
    return (int) bin[bin_idx].size();
}

//-----------------------------------------------------------------------------
void SegmentSearch::query(int bin_idx, std::vector<SegmentNode> &seg_list)
{
    seg_list.insert(seg_list.end(), bin[bin_idx].begin(), bin[bin_idx].end());
}

//-----------------------------------------------------------------------------
void SegmentSearch::query_neighborhood(int bin_idx, 
                                       std::vector<SegmentNode> &seg_list)
{
    query(bin_idx, seg_list);
    for (int j=0;j<bin_nb_neigh[bin_idx];j++) {
        int i = bin_neigh[bin_idx][j];
        query(i, seg_list);
    }
}

//-----------------------------------------------------------------------------
void insertion_sort(int arr[], int n)
{
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

//-----------------------------------------------------------------------------
bool SegmentSearch::is_too_close(double* P, double* Q, int bin_idx, 
                                 double radius)
{
    double radius2 = radius * radius;
    selection.clear();

    // loop over the bins
    for (int j=-1;j<bin_nb_neigh[bin_idx];j++) {
        int i = (j < 0) ? bin_idx : bin_neigh[bin_idx][j];
        // search each bin
        for (size_t k=0;k<bin[i].size();k++) {
            SegmentNode* s = &bin[i][k];
            if (segment_distance2(P, Q, s->P, s->Q) < radius2)
                selection.push_back(s->id);
        }
    }

    // check if the segments in the sphere are connected
    int* rank = selection.data();
    int size = (int) selection.size();
    insertion_sort(rank, size);
    for (int j=1;j<=size;j++)
        if (rank[size-j] != next_id-j)
            return 1;
    return 0;
}

//-----------------------------------------------------------------------------
bool SegmentSearch::too_close_within_a_bin(int idx, double radius)
{
    int n = (int) bin[idx].size();
    SegmentNode* last = &bin[idx][n-1];
    for (int i=0;i<n-1;i++) {
        SegmentNode* s = &bin[idx][i];
        if (segment_distance2(last->P, last->Q, s->P, s->Q) < radius*radius)
            return 1;
    }
    return 0;
}

//-----------------------------------------------------------------------------
void SegmentSearch::find_too_close(double* P, double* Q, int bin_idx, 
                                   double radius, 
                                   std::vector<SegmentNode> &seg_list)
{
    double radius2 = radius * radius;
    seg_list.clear();

    // loop over the bins
    for (int j=-1;j<bin_nb_neigh[bin_idx];j++) {
        int i = (j < 0) ? bin_idx : bin_neigh[bin_idx][j];
        // search each bin
        for (size_t k=0;k<bin[i].size();k++) {
            SegmentNode s = bin[i][k];
            if (segment_distance2(P, Q, s.P, s.Q) < radius2)
                seg_list.push_back(s);
        }
    }

    int size = (int) seg_list.size();
    int *rank = new int [size];
    for (int i=0;i<size;i++)
        rank[i] = seg_list[i].id;
    insertion_sort(rank, size);

    int rank_crit = 0;
    for (int i=size-1;i>=1;i--)
        if (rank[i] - rank[i-1] != 1) {
            rank_crit = rank[i];
            break;
        }

    for (int i=0;i<size;i++)
        if (seg_list[i].id >= rank_crit)
            seg_list[i].id = -1;
    
    delete [] rank;
}
#include <stdlib.h>
#include "segmentdeque.h"
#include "disjointintervals.h"
#include "streamlinetracer.h"

//-----------------------------------------------------------------------------
StreamlineTracer::StreamlineTracer()
{
    nt = 0;
    connect = NULL;
    s_crit = NULL;
    height = NULL;
    admissibles = NULL;
    mask = NULL;
    is_partly_full = NULL;
}

//-----------------------------------------------------------------------------
void StreamlineTracer::initialize(int nt_)
{
    nt = nt_;
    connect = new TriangleConnection [3*nt];
    s_crit = new double [nt];
    height = new double [nt];
    admissibles = new int [nt];
    mask = new DisjointIntervals [nt];
    is_partly_full = new char [nt];
    for (int i=0;i<nt;i++) is_partly_full[i] = 0;
}

//-----------------------------------------------------------------------------
StreamlineTracer::StreamlineTracer(int nt_)
{
    initialize(nt_);
}

//-----------------------------------------------------------------------------
StreamlineTracer::~StreamlineTracer()
{
    delete [] connect;
    delete [] s_crit;
    delete [] height;
    delete [] admissibles;
    delete [] mask;
    delete [] is_partly_full;
}

//-----------------------------------------------------------------------------
void StreamlineTracer::set_connections(int idx, int* idx_adj, 
                                       double* a, double* b)
{
    for (int k=0;k<3;k++) {
        TriangleConnection* c = connect + 3*idx + k;
        c->idx = idx_adj[k];
        c->a = a[k];
        c->b = b[k];
    }
}

//-----------------------------------------------------------------------------
void StreamlineTracer::set_parameters(int idx, double s_crit_, double height_)
{
    s_crit[idx] = s_crit_;
    height[idx] = height_;
}

//-----------------------------------------------------------------------------
void StreamlineTracer::set_region(int size_, int* region_)
{
    region_size = size_;
    region = region_;
}

//-----------------------------------------------------------------------------
int StreamlineTracer::count_admissible_triangles()
{
    nb_admissibles = 0;

    if (region_size > 0) {
        // search triangles near existing streamlines
        // within the specified region
        for (int k=0;k<region_size;k++) {
            int i = region[k];
            admissibles[nb_admissibles] = i;
            nb_admissibles += is_partly_full[i];
        }
        
        // otherwise take any available triangle
        if (nb_admissibles == 0) {
            for (int k=0;k<region_size;k++) {
                int i = region[k];
                admissibles[nb_admissibles] = i;
                nb_admissibles += !mask[i].is_full();
            }
        }
    } else {
        // search triangles near existing streamlines
        for (int i=3;i<nt;i+=4) {
            admissibles[nb_admissibles] = i-3;
            nb_admissibles += is_partly_full[i-3];
            admissibles[nb_admissibles] = i-2;
            nb_admissibles += is_partly_full[i-2];
            admissibles[nb_admissibles] = i-1;
            nb_admissibles += is_partly_full[i-1];
            admissibles[nb_admissibles] = i;
            nb_admissibles += is_partly_full[i];
        }
        for (int i=(nt>>2)<<2;i<nt;i++) {
            admissibles[nb_admissibles] = i;
            nb_admissibles += is_partly_full[i];
        }

        // otherwise take any available triangle
        if (nb_admissibles == 0) {
            for (int i=0;i<nt;i++) {
                admissibles[nb_admissibles] = i;
                nb_admissibles += !mask[i].is_full();
            }
        }
    }

    // if there is no more space, return 0 (fail)
    return nb_admissibles;
}

//-----------------------------------------------------------------------------
int StreamlineTracer::pick_random(int n, int* idx, double* s, int* rnd_stream)
{
    for (int i=0;i<n;i++) {
        idx[i] = admissibles[rnd_stream[2*i] % nb_admissibles];
        s[i] = mask[idx[i]].pick_random(rnd_stream[2*i+1]);
    }
    return n;
}

//-----------------------------------------------------------------------------
double StreamlineTracer::trace_streamline(int idx_seed, double s_seed, 
                                          SegmentDeque* queue)
{
    TriangleConnection* c;
    int prev_idx, idx = idx_seed;
    double s = s_seed;
    double sc, L = 0;

    #define UPDATE_STEP \
        prev_idx = idx; \
        idx = c->idx; \
        s = c->a * s + c->b;
    
    #define UPDATE_LENGTH \
        sc = s_crit[idx]; \
        L += height[idx] * ((s <= sc) ? s/sc : (1-s)/(1-sc));
    
    #define FIRST_STEP_FORWARD \
        c = connect + 3*idx + 1 + (s > s_crit[idx]); \
        UPDATE_STEP;

    #define FIRST_STEP_BACKWARD \
        c = connect + 3*idx; \
        UPDATE_STEP;

    #define NEXT_STEP \
        c = connect + 3*idx; \
        c += (c->idx == prev_idx) * (1 + (s > s_crit[idx])); \
        UPDATE_STEP;
    
    queue->clear();
    queue->append(idx, s);
    UPDATE_LENGTH;

    // forward integration
    FIRST_STEP_FORWARD;
    while (idx != -1) {
        UPDATE_LENGTH;
        if (!queue->append(idx, s)) break;
        if (mask[idx].contains(s)) break;
        NEXT_STEP;
    }
    // backward integration
    queue->peek_first(idx, s);
    FIRST_STEP_BACKWARD;
    while (idx != -1) {
        UPDATE_LENGTH;
        if (!queue->prepend(idx, s)) break;
        if (mask[idx].contains(s)) break;
        NEXT_STEP;
    }
    return L;

    #undef UPDATE_STEP
    #undef UPDATE_LENGTH
    #undef NEXT_STEP
    #undef FIRST_STEP_FORWARD
    #undef FIRST_STEP_BACKWARD
}

//-----------------------------------------------------------------------------
inline void StreamlineTracer::set_mask(int idx, double s1, double s2)
{
    mask[idx].insert(s1, s2);
    bool state = !mask[idx].is_full() && !mask[idx].is_empty();
    is_partly_full[idx] = state;
}
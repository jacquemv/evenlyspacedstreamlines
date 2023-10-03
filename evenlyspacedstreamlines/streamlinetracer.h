#ifndef STREAMLINETRACER_H_
#define STREAMLINETRACER_H_

//-----------------------------------------------------------------------------
typedef struct TriangleConnection {
    int idx;
    double a, b; // s_next = a * s_prev + b
} TriangleConnection;

//-----------------------------------------------------------------------------
class StreamlineTracer {
public:
    int nt;
    TriangleConnection* connect; // 3*nt
    double *s_crit, *height;
    DisjointIntervals *mask;

    StreamlineTracer();
    StreamlineTracer(int nt_);
    ~StreamlineTracer();
    void initialize(int nt_);

    int count_admissible_triangles();
    int pick_random(int n, int* idx, double* s, int* rnd_stream);
    double trace_streamline(int idx, double s, SegmentDeque* queue);
    inline void set_mask(int idx, double s1, double s2);

    void set_connections(int idx, int* idx_adj, double* a, double* b);
    void set_parameters(int idx, double s_crit_, double height_);
    void set_region(int size, int* region);

//private:
    int nb_admissibles;
    int* admissibles;
    char* is_partly_full;
    int region_size;
    int* region;
};

#endif
#ifndef ENGINE_H_
#define ENGINE_H_

#include "algebra.h"
#include "curvature.h"
#include "disjointintervals.h"
#include "geometry.h"
#include "quadraticforms.h"
#include "segmentdeque.h"
#include "segmentsearch.h"
#include "streamlinecontainer.h"
#include "streamlinetracer.h"
#include "triangle.h"
#include "triangularmesh.h"
#include <vector>


//-----------------------------------------------------------------------------
class StreamlineEngine {
public:
    double radius, max_angle;
    int max_nb_seeds;
    int avoid_u_turns;
    int parallel;
    bool oriented_streamlines;

    StreamlineEngine();
    ~StreamlineEngine();

    void initialize(int nv_, int nt_, double* ver_, int* tri_, 
                    double* orient_, int orthogonal, 
                    int allow_tweaking_orientation);
    // max_angle = maximum angle between successive segment of a streamline
    //             default = 90 degrees ; > 90 = relieve the constraint ; 
    //             < 90 = very strict constraint ; 180 = deactivated
    void setup(double radius_, int max_length, int max_nb_seeds_, 
               int avoid_u_turns, double max_angle, 
               int oriented_streamlines_,
               double singularity_mask_radius,
               unsigned int random_seed, int parallel_, int num_threads);
    // call after setup()
    void define_seed_region(int seed_region_size, int* seed_region);
    void run();

    double min_altitude, max_base;
    int euler_characteristic;
    double neighborhood_mean_size;
    unsigned int saved_random_seed;
    int error_code;

    int generate_seed_points();
    SegmentDeque* optimize_streamline(int nb_seeds);
    SegmentDeque* optimize_streamline_parallel(int nb_seeds);
    int cut_nearly_periodic(SegmentDeque* queue);
    int cut_nearly_periodic_extremity(SegmentDeque* queue);
    int cut_highly_curved(SegmentDeque* queue);
    int cut_extremities(SegmentDeque* queue);
    void update_mask(SegmentDeque* queue);
    void update_mask_parallel(SegmentDeque* queue);
    void mask_sphere(int idver, double R);
    void store_streamline(SegmentDeque* queue, 
                          std::vector<StreamlineContainer*> &storage);
    int get_nb_streamlines();
    int get_streamline_size(int no);
    void get_streamline(int no, double* points, int* idx);
    void get_lengths(double* L);
    
    void print_timings();

//private:
    TriangularMesh mesh;
    StreamlineTracer tracer;
    std::vector<StreamlineContainer*> collection;
    SegmentSearch singlecurve, allsegments;
    std::vector<SegmentNode> segment_list;
    Curvature curvature;

    void setup_geometry(double max_angle);
    void setup_masking_singularities(double radius_ratio);
    void setup_masking_invalid_orientation();
    void setup_streamline_storage(int maxlen);
    void setup_random_number_generator(int random_seed);

    #ifdef DEBUG_STREAMLINES
    std::vector<StreamlineContainer*> collection_step1, collection_step2, 
                                      collection_step3, collection_step4;
    void store_uncut_streamline(SegmentDeque* queue, 
                                std::vector<StreamlineContainer*> &storage);
    void print_streamline(SegmentDeque* queue);
    void check_streamline_size(const char * mesg);
    char* mask_filename;
    void write_mask_polygons();
    void write_orientation_vectors(const char* fname, int n);
    void write_lines_of_block(const char* fname);
    #endif

private:
    int max_threads;
    SegmentDeque **current_streamline, **thread_best_streamline, 
                  *best_streamline;
    double* best_length;
    int* idx_seed;
    double* s_seed;
    int* rnd_stream;
};

#endif
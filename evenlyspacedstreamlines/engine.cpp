#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cassert>
#include <omp.h>

#include <vector>
#include <fstream>
#include <iostream>
#include "engine.h"

#include "algebra.cpp"
#include "curvature.cpp"
#include "disjointintervals.cpp"
#include "geometry.cpp"
#include "quadraticforms.cpp"
#include "segmentdeque.cpp"
#include "segmentsearch.cpp"
#include "streamlinecontainer.cpp"
#include "streamlinetracer.cpp"
#include "triangle.cpp"
#include "triangularmesh.cpp"

#define CACHE_LINE 8 // * sizeof(double) bytes
#define MAXLENGTH_MINVAL 256

#ifdef TIMING
    #include "timing.h"
    #define TIMING_ON(n) timing_on(n)
    #define TIMING_OFF(n) timing_off(n)
#else
    #define TIMING_ON(n)  //printf("begin %d\n", n);
    #define TIMING_OFF(n) //printf("end %d\n", n);
#endif

//-----------------------------------------------------------------------------
StreamlineEngine::StreamlineEngine()
{
    max_threads = 0;
    current_streamline = NULL;
    thread_best_streamline = NULL, 
    best_streamline = NULL;
    best_length = NULL;
    idx_seed = NULL;
    s_seed = NULL;
    rnd_stream = NULL;
}

//-----------------------------------------------------------------------------
void StreamlineEngine::initialize(int nv_, int nt_, double* ver_, int* tri_, 
                                   double* orient_, int orthogonal,
                                   int allow_tweaking_orientation)
{
    TIMING_ON(31); TIMING_ON(0);
    error_code = 0;
    int err = mesh.initialize(nv_, nt_, ver_, tri_, orient_, orthogonal,
                              allow_tweaking_orientation);
    if (err < 0) {
        error_code = 3;
    }
    #ifdef DEBUG_STREAMLINES
    mask_filename = (char*) NULL;
    #endif
    TIMING_OFF(0); TIMING_OFF(31);
}

//-----------------------------------------------------------------------------
StreamlineEngine::~StreamlineEngine()
{
    for (int i=0;i<max_threads;i++) {
        delete current_streamline[CACHE_LINE*i];
        delete thread_best_streamline[CACHE_LINE*i];
    }
    delete [] current_streamline;
    delete [] thread_best_streamline;
    delete [] best_length;
    delete [] idx_seed;
    delete [] s_seed;
    delete [] rnd_stream;
    for (unsigned int i=0;i<collection.size();i++) 
        delete collection[i];
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup(double radius_, int max_length, int max_nb_seeds_,
                             int avoid_u_turns_, double max_angle, 
                             int oriented_streamlines_,
                             double singularity_mask_radius,
                             unsigned int random_seed, 
                             int parallel_, int num_threads)
{
    TIMING_ON(31); TIMING_ON(0);
    radius = radius_;
    max_nb_seeds = max_nb_seeds_;
    avoid_u_turns = avoid_u_turns_;
    oriented_streamlines = oriented_streamlines_;

    parallel = parallel_;
    max_threads = parallel ? omp_get_max_threads() : 1;
    if (num_threads > 0) {
        if (num_threads > max_threads) num_threads = max_threads;
        omp_set_num_threads(num_threads);
    }
    
    int default_max_length = (mesh.nt < MAXLENGTH_MINVAL) ? MAXLENGTH_MINVAL 
                                                           : mesh.nt;
    int streamline_max_length = (max_length > 0) ? max_length 
                                                 : default_max_length; 
    if (avoid_u_turns) {
        if (max_angle > 90) max_angle = 90;
    }
    setup_geometry(max_angle);
    setup_masking_singularities(singularity_mask_radius);
    setup_masking_invalid_orientation();
    setup_streamline_storage(streamline_max_length);
    setup_random_number_generator(random_seed);

    singlecurve.initialize(mesh.nt, mesh.t_nb_extneigh_tri, 
                           mesh.t_extneigh_tri);
    allsegments.initialize(mesh.nt, mesh.t_nb_extneigh_tri, 
                           mesh.t_extneigh_tri);
    curvature.initialize(streamline_max_length);
    TIMING_OFF(0); TIMING_OFF(31);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::define_seed_region(int seed_region_size, 
                                          int* seed_region)
{
    tracer.set_region(seed_region_size, seed_region);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup_geometry(double max_angle)
{
    mesh.compute_neighborhood(radius, parallel);
    double cos_limit_angle = cos(max_angle*0.017453292519943295); // pi/180

    tracer.initialize(mesh.nt);
    tracer.set_region(0, NULL);
    for (int i=0;i<mesh.nt;i++) {
        int adj[3];
        double a[3], b[3], cos_angle[3];
        bool is_continuous[3];
        mesh.get_triangle_connections(i, adj, a, b, is_continuous, cos_angle);
        for (int k=0;k<3;k++) {
            if (cos_angle[k] < cos_limit_angle)
                adj[k] = -1;
        }
        if (oriented_streamlines) {
            for (int k=0;k<3;k++) {
                if (!is_continuous[k])
                    adj[k] = -1;
            }
        }
        tracer.set_connections(i, adj, a, b);
        double s_crit, height;
        mesh.get_triangle_parameters(i, s_crit, height);
        tracer.set_parameters(i, s_crit, height);
        // orientation parallel to an edge
        if (s_crit == 0 || s_crit == 1) error_code = 2;
    }

    min_altitude = mesh.triangle_min_altitude;
    max_base = mesh.triangle_max_base;
    euler_characteristic = mesh.euler_characteristic;
    neighborhood_mean_size = mesh.neighborhood_mean_size;

    if (min_altitude <= 0) // triangle area = 0
        error_code = 1;
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup_masking_singularities(double radius_ratio)
{
    if (radius_ratio <= 0) return;

    std::vector<int> focus, uturn, node;
    mesh.find_singularities(focus, uturn, node);

    for (size_t i=0;i<focus.size();i++)
        mask_sphere(focus[i], radius_ratio * radius);

    for (size_t i=0;i<uturn.size();i++)
        mask_sphere(uturn[i], radius_ratio * radius);

    for (size_t i=0;i<node.size();i++)
        mask_sphere(node[i], radius_ratio * radius);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup_masking_invalid_orientation()
{
    for (int i=0;i<mesh.nt;i++) {
        double *v = mesh.orient+3*i;
        if (((v[0] == 0) + (v[1] == 0) + (v[2] == 0) == 3) ||
            (isfinite(v[0]) + isfinite(v[1]) + isfinite(v[2]) != 3)) {
                tracer.set_mask(i, 0, 1);
            }
    }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup_streamline_storage(int maxlen)
{
    int nblock = parallel ? CACHE_LINE * max_threads : 1;

    current_streamline = new SegmentDeque* [nblock];
    thread_best_streamline = new SegmentDeque* [nblock];
    best_length = new double [nblock];

    for (int i=0;i<max_threads;i++) {
        current_streamline[CACHE_LINE*i] = new SegmentDeque(maxlen);
        thread_best_streamline[CACHE_LINE*i] = new SegmentDeque(maxlen);
    }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::setup_random_number_generator(int random_seed)
{
    if (random_seed == 0)
        saved_random_seed = (unsigned int) time(NULL);
    else
        saved_random_seed = random_seed;
    srand(saved_random_seed);

    rnd_stream = new int [2*max_nb_seeds];
    idx_seed = new int [max_nb_seeds];
    s_seed = new double [max_nb_seeds];
}

//-----------------------------------------------------------------------------
#ifdef DEBUG_STREAMLINES
    #define DEBUG_RUN_STEP0 
    #define DEBUG_RUN_STEP1 store_uncut_streamline(queue, collection_step1);
    #define DEBUG_RUN_STEP2 store_uncut_streamline(queue, collection_step2);
    #define DEBUG_RUN_STEP3 store_uncut_streamline(queue, collection_step3);
    #define DEBUG_RUN_STEP4 store_streamline(queue, collection_step4);
    #define DEBUG_RUN_STEP5 write_mask_polygons();
    #define DEBUG_RUN_STEP6
#else
    #define DEBUG_RUN_STEP0
    #define DEBUG_RUN_STEP1
    #define DEBUG_RUN_STEP2
    #define DEBUG_RUN_STEP3
    #define DEBUG_RUN_STEP4
    #define DEBUG_RUN_STEP5
    #define DEBUG_RUN_STEP6
#endif

//-----------------------------------------------------------------------------
void StreamlineEngine::run()
{
    DEBUG_RUN_STEP0;
    TIMING_ON(31); 
    SegmentDeque* queue;
    while (1) {
        TIMING_ON(1);
        int nb_seeds = generate_seed_points();
        if (!nb_seeds) break;
        TIMING_OFF(1);
        
        TIMING_ON(2);
        if (parallel)
            queue = optimize_streamline_parallel(nb_seeds);
        else
            queue = optimize_streamline(nb_seeds);
        if (!queue) break;
        queue->set_data(0, 1, 0, 1);
        TIMING_OFF(2);
        DEBUG_RUN_STEP1;

        TIMING_ON(3);
        if (cut_nearly_periodic(queue)) {
            cut_nearly_periodic_extremity(queue);
        }
        TIMING_OFF(3);
        DEBUG_RUN_STEP2;

        TIMING_ON(4);
        if (avoid_u_turns) 
            cut_highly_curved(queue);
        TIMING_OFF(4);
        DEBUG_RUN_STEP3;

        TIMING_ON(5);
            cut_extremities(queue);
        TIMING_OFF(5);
        DEBUG_RUN_STEP4;

        TIMING_ON(6);
        if (parallel)
            update_mask_parallel(queue);
        else
            update_mask(queue);
        TIMING_OFF(6);
        DEBUG_RUN_STEP5;

        TIMING_ON(7);
        store_streamline(queue, collection);
        TIMING_OFF(7);
        DEBUG_RUN_STEP6;
    }
    TIMING_OFF(31); 
}

//-----------------------------------------------------------------------------
void StreamlineEngine::print_timings()
{
    #ifdef TIMING
    timing_print(0, "setup");
    timing_print(1, "seed points");
    timing_print(10, "  admissible");
    timing_print(11, "  random");
    timing_print(12, "  pick");
    timing_print(2, "streamlines");
    timing_print(3, "nearly-periodic");
    timing_print(4, "highly-curved");
    timing_print(5, "extremities");
    timing_print(6, "mask");
    timing_print(7, "storage");
    timing_print(31, "total");
    #endif
}

//-----------------------------------------------------------------------------
int StreamlineEngine::generate_seed_points()
{
    TIMING_ON(10);
    int count = tracer.count_admissible_triangles();
    TIMING_OFF(10);
    if (count == 0) return 0;

    TIMING_ON(11);
    for (int i=0;i<2*max_nb_seeds;i++)
        rnd_stream[i] = rand();
    TIMING_OFF(11);
    
    TIMING_ON(12);
    int nb_seeds = tracer.pick_random(max_nb_seeds, idx_seed, s_seed, 
                                      rnd_stream);
    TIMING_OFF(12);

    return nb_seeds;
}

//-----------------------------------------------------------------------------
#define SWAP_STREAMLINES(S1, S2) { \
    SegmentDeque* temp = S1; S1 = S2; S2 = temp; \
}

//-----------------------------------------------------------------------------
SegmentDeque* StreamlineEngine::optimize_streamline(int nb_seeds)
{
    double L, Lmax = 0;
    for (int i=0;i<nb_seeds;i++) {
        L = tracer.trace_streamline(idx_seed[i], s_seed[i], 
                                    current_streamline[0]);
        if (L > Lmax) {
            Lmax = L;
            SWAP_STREAMLINES(current_streamline[0], thread_best_streamline[0]);
        }
    }
    best_streamline = thread_best_streamline[0];
    return best_streamline;
}

//-----------------------------------------------------------------------------
SegmentDeque* StreamlineEngine::optimize_streamline_parallel(int nb_seeds)
{
    int nthreads;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp master
        nthreads = omp_get_num_threads();


        int lane = CACHE_LINE * thread_id;
        SegmentDeque* priv_current_streamline = current_streamline[lane];
        SegmentDeque* priv_best_streamline = thread_best_streamline[lane];
        double L, Lmax = 0;

        #pragma omp for schedule(static,16)
        for (int i=0;i<nb_seeds;i++) {
            L = tracer.trace_streamline(idx_seed[i], s_seed[i],
                                        priv_current_streamline);
            if (L > Lmax) {
                Lmax = L;
                SWAP_STREAMLINES(priv_current_streamline, 
                                 priv_best_streamline);
            }
        }

        best_length[lane] = Lmax;
        thread_best_streamline[lane] = priv_best_streamline;
        current_streamline[lane] = priv_current_streamline;
    }

    // combine threads
    double Lmax = best_length[0];
    best_streamline = thread_best_streamline[0];
    for (int i=1;i<nthreads;i++) {
        int lane = CACHE_LINE*i;
        if (best_length[lane] > Lmax) {
            Lmax = best_length[lane];
            best_streamline = thread_best_streamline[lane];
        }
    }
    return best_streamline;
}

//-----------------------------------------------------------------------------
#define QUEUE_FOR_EACH_ELEMENT(queue, idx, s) \
    queue->rewind(); \
    int idx=0; double s=0; \
    while (queue->peek_next(idx, s))

#define GET_SEGMENT(mesh, idx, s, P, Q) \
    double P[3], Q[3]; \
    mesh.get_segment(idx, s, P, Q);

#define GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q) \
    double P[3], Q[3]; \
    if (queue->is_extremity()) { \
        double t1, t2; \
        queue->get_data(t1, t2); \
        mesh.get_partial_segment(idx, s, t1, t2, P, Q); \
    } else mesh.get_segment(idx, s, P, Q);

#define GET_PARTIAL_SEGMENT_PARALLEL(mesh, queue, idx, s, P, Q) \
    double P[3], Q[3]; \
    if (pos == queue->first || pos == queue->last) { \
        double t1, t2; \
        queue->get_data(pos, t1, t2); \
        mesh.get_partial_segment(idx, s, t1, t2, P, Q); \
    } else mesh.get_segment(idx, s, P, Q);

#define GET_SEGMENT_NO(mesh, queue, seg_no, idx, s, P, Q) \
    double P[3], Q[3]; \
    queue->peek_position(seg_no, idx, s); \
    mesh.get_segment(idx, s, P, Q);

//-----------------------------------------------------------------------------
int StreamlineEngine::cut_highly_curved(SegmentDeque* queue)
{
    curvature.reset();
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_SEGMENT(mesh, idx, s, P, Q);
        curvature.append_segment(P, Q);
    }
    int start, len;
    curvature.low_curvature_interval(radius, start, len);
    int cut = queue->size() - len;
    queue->select(start, len);
    return cut;
}

//-----------------------------------------------------------------------------
int StreamlineEngine::cut_nearly_periodic(SegmentDeque* queue)
{
    int len = 0;
    singlecurve.clear();
    if (queue->prepend_failed() && !queue->append_failed())
        queue->flip();
    int stop = 0;
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_SEGMENT(mesh, idx, s, P, Q);
        stop = singlecurve.is_too_close(P, Q, idx, radius);
        len++;
        if (stop) break;
        int binsize = singlecurve.insert(P, Q, idx);
        if (binsize  > 1) {
            // avoid long loops when the entire streamline is in a sphere 
            // of radius r
            if (singlecurve.too_close_within_a_bin(idx, radius)) {
                len--;
                break;
            }
        }
    }
    queue->select(0, len);
    return stop;
}

//-----------------------------------------------------------------------------
int StreamlineEngine::cut_nearly_periodic_extremity(SegmentDeque* queue)
{
    if (queue->size() < 2)
        return 0;

    int last_no = -1, before_last_no = -2; // cut only the 'last' end
    int idx = 0; double s = 0;
    GET_SEGMENT_NO(mesh, queue, before_last_no, idx, s, P2, Q2);
    GET_SEGMENT_NO(mesh, queue, last_no, idx, s, P, Q);
    int code = connect_segments(P, Q, P2, Q2);

    double tmin = 1, tmax = 0, t1, t2;
    #define UPDATE_INTERVAL { if (t1<tmin) tmin = t1; \
                              if (t2>tmax) tmax = t2; }
    std::vector<SegmentNode> seg_list;
    singlecurve.find_too_close(P, Q, idx, radius, seg_list);

    for (unsigned int i=0;i<seg_list.size();i++) {
        SegmentNode *s = &seg_list[i];
        if (s->id < 0) continue; // segments tagged by find_too_close
        if (intersect_sphere(P, Q, s->P, radius, t1, t2))
            UPDATE_INTERVAL;
        if (intersect_sphere(P, Q, s->Q, radius, t1, t2))
            UPDATE_INTERVAL;
        if (intersect_cylinder(P, Q, s->P, s->Q, radius, t1, t2))
            UPDATE_INTERVAL;
    }

    if (code / 2) { // if P is the last point of the streamline
        t1 = tmax; t2 = 1; if (t1 > 1) t1 = 1;
    } else { // Q is the last point
        t1 = 0; t2 = tmin; if (t2 < 0) t2 = 0;
    }
    queue->data_last[0] = t1;
    queue->data_last[1] = t2;
    return 1;
}

//-----------------------------------------------------------------------------
int StreamlineEngine::cut_extremities(SegmentDeque* queue)
{
    if (queue->size() < 2)
        return 0;

    int last_no = 0, before_last_no = 1;
    for (int end_point=0;end_point<2;end_point++) {

        int idx = 0; double s = 0;
        GET_SEGMENT_NO(mesh, queue, before_last_no, idx, s, P2, Q2);
        GET_SEGMENT_NO(mesh, queue, last_no, idx, s, P, Q);
        int code = connect_segments(P, Q, P2, Q2);

        double tmin = 1, tmax = 0, t1, t2;
        #define UPDATE_INTERVAL { if (t1<tmin) tmin = t1; \
                                  if (t2>tmax) tmax = t2; }
        segment_list.clear();
        allsegments.query_neighborhood(idx, segment_list);
        for (unsigned int i=0;i<segment_list.size();i++) {
            SegmentNode *s = &segment_list[i];
            if (intersect_sphere(P, Q, s->P, radius, t1, t2))
                UPDATE_INTERVAL;
            if (intersect_sphere(P, Q, s->Q, radius, t1, t2))
                UPDATE_INTERVAL;
            if (intersect_cylinder(P, Q, s->P, s->Q, radius, t1, t2))
                UPDATE_INTERVAL;
        }

        if (code / 2) { // if P is the last point of the streamline
            t1 = tmax; t2 = 1; if (t1 > 1) t1 = 1;
        } else { // Q is the last point
            t1 = 0; t2 = tmin; if (t2 < 0) t2 = 0;
        }
        if (end_point) {
            if (t1 > queue->data_last[0]) // compute the intersection with
                queue->data_last[0] = t1; // the current interval
            if (t2 < queue->data_last[1]) // e.g. from cut_nearly_periodic
                queue->data_last[1] = t2;
        } else{
            if (t1 > queue->data_first[0])
                queue->data_first[0] = t1;
            if (t2 < queue->data_first[1])
                queue->data_first[1] = t2;
        }

        last_no = -1, before_last_no = -2;
    }
    return 0;
}

//-----------------------------------------------------------------------------
void StreamlineEngine::update_mask(SegmentDeque* queue)
{
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        // mask of the same triangle
        double ds = radius * mesh.shrink_factor(idx);
        tracer.set_mask(idx, s - ds, s + ds);

        // mask of the neighboring triangles
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        int nb_neigh, *neigh_list;
        mesh.get_neighborhood(idx, neigh_list, nb_neigh);
        for (int k=0;k<nb_neigh;k++) {
            int idx_neigh = neigh_list[k];
            double s1, s2;
            mesh.compute_mask(idx_neigh, P, Q, radius, s1, s2);
            tracer.set_mask(idx_neigh, s1, s2);
        }

        // update the list fo all segments
        allsegments.insert(P, Q, idx);
    }
}

//-----------------------------------------------------------------------------
typedef struct MaskInstruction {
    int idx;
    double s1, s2;
} MaskInstruction;

void StreamlineEngine::update_mask_parallel(SegmentDeque* queue)
{
    std::vector<MaskInstruction> *instructions;
    int nthreads, max_threads = omp_get_max_threads();
    instructions = new std::vector<MaskInstruction> [max_threads*CACHE_LINE];

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int lane = thread_id * CACHE_LINE;
        #pragma omp master
        nthreads = omp_get_num_threads();
        #pragma omp for
        for (int pos=queue->first; pos<=queue->last; pos++) {
            int idx = queue->array[pos].idx;
            double s = queue->array[pos].s;

            // mask of the same triangle
            double ds = radius * mesh.shrink_factor(idx);
            MaskInstruction instr = {idx, s-ds, s+ds};
            instructions[lane].push_back(instr);

            // mask of the neighboring triangles
            GET_PARTIAL_SEGMENT_PARALLEL(mesh, queue, idx, s, P, Q);
            int nb_neigh, *neigh_list;
            mesh.get_neighborhood(idx, neigh_list, nb_neigh);
            for (int k=0;k<nb_neigh;k++) {
                int idx_neigh = neigh_list[k];
                double s1, s2;
                mesh.compute_mask(idx_neigh, P, Q, radius, s1, s2);
                MaskInstruction instr = {idx_neigh, s1, s2};
                instructions[lane].push_back(instr);
            }
        }
    } // end parallel

    // execute instructions
    for (int k=0;k<nthreads;k++) {
        int lane = k * CACHE_LINE;
        for (unsigned int i=0;i<instructions[lane].size();i++) {
            tracer.set_mask(instructions[lane][i].idx, 
                            instructions[lane][i].s1, 
                            instructions[lane][i].s2);
        }
    }
    delete [] instructions;

    // update the list fo all segments
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        allsegments.insert(P, Q, idx);
    }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::mask_sphere(int idver, double R)
{
    int nb_neigh, *neigh_list;
    double *C = mesh.ver + 3*idver;
    mesh.get_incident_triangles(idver, neigh_list, nb_neigh);
    for (int k=0;k<nb_neigh;k++) {
        int idx_neigh = neigh_list[k];
        double s1, s2;
        mesh.face[idx_neigh].coordinate_mask_sphere(C, R, s1, s2);
        tracer.set_mask(idx_neigh, s1, s2);
    }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::store_streamline(SegmentDeque* queue, 
                                std::vector<StreamlineContainer*> &storage)
{
    int nb_points = queue->size() + 1;
    StreamlineContainer *streamline = new StreamlineContainer(nb_points);
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        streamline->push_segment(P, Q, idx);
    }
    if (oriented_streamlines) {
        double vector[3];
        int id_tri = streamline->first_vector(vector);
        if ((vdot(vector, mesh.face[id_tri].orient) < 0) 
                ^ mesh.face[id_tri].reverse_orient)
            streamline->reverse();
    }
    storage.push_back(streamline);
}

//-----------------------------------------------------------------------------
int StreamlineEngine::get_nb_streamlines()
{
    return (int) collection.size();
}

//-----------------------------------------------------------------------------
int StreamlineEngine::get_streamline_size(int no)
{
    return collection[no]->n;
}

//-----------------------------------------------------------------------------
void StreamlineEngine::get_streamline(int no, double* points, int* idx)
{
    int n = collection[no]->n;
    for (int i=0;i<3*n;i++)
        points[i] = collection[no]->points[i];
    for (int i=0;i<n-1;i++)
        idx[i] = collection[no]->idx[i];
}

//-----------------------------------------------------------------------------
void StreamlineEngine::get_lengths(double* L)
{
    for (unsigned int i=0;i<collection.size();i++)
        L[i] = collection[i]->length();
}

//=============================================================================
#ifdef DEBUG_STREAMLINES
//-----------------------------------------------------------------------------
void StreamlineEngine::store_uncut_streamline(SegmentDeque* queue, 
                                std::vector<StreamlineContainer*> &storage)
{
    int nb_points = queue->size() + 1;
    StreamlineContainer *streamline = new StreamlineContainer(nb_points);
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s) {
        GET_SEGMENT(mesh, idx, s, P, Q);
        streamline->push_segment(P, Q, idx);
    }
    storage.push_back(streamline);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::print_streamline(SegmentDeque* queue)
{
    QUEUE_FOR_EACH_ELEMENT(queue, idx, s){
        GET_PARTIAL_SEGMENT(mesh, queue, idx, s, P, Q);
        printf("%.7f %.7f %.7f   %.7f %.7f %.7f\n", 
               P[0], P[1], P[2], Q[0], Q[1], Q[2]);
    }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::check_streamline_size(const char * mesg)
{
    for (unsigned int i=0;i<collection.size();i++)
        if (collection[i]->n < 0) {
            printf("%s: %d %d\n", mesg, i, collection[i]->n);
            exit(1);
        }
}

//-----------------------------------------------------------------------------
void StreamlineEngine::write_mask_polygons()
{
    char buffer[256], default_filename[32] = "mask";
    if (mask_filename)
        sprintf(buffer, "%s%d.txt", mask_filename, (int)collection.size());
    else
        sprintf(buffer, "%s%d.txt", default_filename, (int)collection.size());
    FILE* fid = fopen(buffer, "w+");
    #define PUSH(A) fprintf(fid, "%.5f %.5f %.5f ", A[0], A[1], A[2]);
    for (int idx=0;idx<mesh.nt;idx++) {
        for (int i=0;i<tracer.mask[idx].n;i++) {
            double s1 = tracer.mask[idx].t[2*i];
            double s2 = tracer.mask[idx].t[2*i+1];
            double P[3], Q[3];
            mesh.get_segment(idx, s1, P, Q); PUSH(P); PUSH(Q);
            double sc = tracer.s_crit[idx];
            if (s1 < sc && sc < s2) {
                mesh.get_segment(idx, sc, P, Q); PUSH(Q);
            }
            mesh.get_segment(idx, s2, P, Q); PUSH(Q); PUSH(P);
        }
        if (tracer.mask[idx].n > 0) fprintf(fid, "\n");
    }
    #undef PUSH
    fclose(fid);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::write_orientation_vectors(const char* fname, int n)
{
    FILE* fid = fopen(fname, "w+");
    for (int i=0;i<mesh.nt;i++) {
        for (int j=0;j<n;j++) {
            double P[6];
            mesh.get_segment(i, ((double)j+1.0)/(n+1.0), P, P+3);
            for (int k=0;k<6;k++) fprintf(fid, "%.5f ", P[k]);
            fprintf(fid, "\n");
        }
    }
    fclose(fid);
}

//-----------------------------------------------------------------------------
void StreamlineEngine::write_lines_of_block(const char* fname)
{
    FILE* fid = fopen(fname, "w+");
    for (int i=0;i<mesh.nt;i++) {
        int i1, i2;
        if (tracer.connect[3*i].idx == -1) {
            i1 = mesh.face[i].idx_origin;
            i2 = mesh.face[i].idx_base;
            fprintf(fid, "%d %d\n", i1, i2);
        }
        if (tracer.connect[3*i+1].idx == -1) {
            i1 = mesh.face[i].idx_origin;
            i2 = mesh.face[i].idx_apex;
            fprintf(fid, "%d %d\n", i1, i2);
        }
        if (tracer.connect[3*i+2].idx == -1) {
            i1 = mesh.face[i].idx_apex;
            i2 = mesh.face[i].idx_base;
            fprintf(fid, "%d %d\n", i1, i2);
        }
    }
    fclose(fid);
}

#endif
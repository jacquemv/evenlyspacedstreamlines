#include "engine.cpp"
#include "readmesh.cpp"
#include "tictoc.h"

int main()
{
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off("../data/test_finemesh2d.off", nv, nt, ver, tri);
    read_orient_txt("../data/test_finemesh2d.orient", nt, orient);

    StreamlineEngine engine;
    engine.initialize(nv, nt, ver, tri, orient, 0, 0);
    double radius = 0.01;

    /* setup arguments:
    double radius, int max_length, int max_nb_seeds, 
    int avoid_u_turns, double max_angle, 
    double singularity_mask_radius,
    unsigned int random_seed, int parallel, int num_threads */
    engine.setup(radius, nt, 32, 1, 90, 0, 0.1, 129873, 1, -1);
    engine.run();
    fprintf(stderr, "%d streamlines\n", engine.get_nb_streamlines());
    fprintf(stderr, "radius/altitude = %.2f\n", radius/engine.min_altitude);
    engine.print_timings();
    
    for (int i=0;i<engine.collection.size();i++) {
        engine.collection[i]->print();
        printf("\n");
    }
}
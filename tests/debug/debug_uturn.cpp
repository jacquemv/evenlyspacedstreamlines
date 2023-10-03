#define DEBUG_STREAMLINES
#include "engine.cpp"
#include "readmesh.cpp"
#include "debug.cpp"


int main()
{
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off("../data/test_uturn.off", nv, nt, ver, tri);
    read_orient_txt("../data/test_uturn.orient", nt, orient);

    StreamlineEngine engine;
    engine.initialize(nv, nt, ver, tri, orient, 0, 1);
    double radius = 0.2;

    /* setup arguments:
    double radius, int max_length, int max_nb_seeds, 
    int avoid_u_turns, double max_angle, 
    int oriented_streamlines,
    double singularity_mask_radius,
    unsigned int random_seed, int parallel */
    engine.setup(radius, 0, 32, 
                 1, 90, 
                 0, 
                 -0.1, 
                 1656514220UL, 0);

    char fname[] = "output/mask";
    engine.mask_filename = fname;
    engine.write_orientation_vectors("output/orient.txt", 4);
    engine.write_lines_of_block("output/blocks.txt");
    engine.run();

    write_streamline_collection("output/streamlines1.txt", engine.collection_step1);
    write_streamline_collection("output/streamlines2.txt", engine.collection_step2);
    write_streamline_collection("output/streamlines3.txt", engine.collection_step3);
    write_streamline_collection("output/streamlines4.txt", engine.collection_step4);
    write_streamline_collection("output/streamlines5.txt", engine.collection);
}
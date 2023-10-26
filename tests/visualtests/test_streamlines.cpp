#include "engine.cpp"
#include "readmesh.cpp"

int main()
{
    int nv, nt;
    double *ver, *orient;
    int *tri;
    char mesh_file[] = "../data/test_mesh2d.off";
    char orient_file[] = "../data/test_mesh2d.orient2";
    //char mesh_file[] = "../data/test_annulus.off";
    //char orient_file[] = "../data/test_annulus.orient";

    read_off(mesh_file, nv, nt, ver, tri);
    read_orient_txt(orient_file, nt, orient);

    StreamlineEngine engine;
    engine.initialize(nv, nt, ver, tri, orient, 0, 0);
    double radius = 0.1;
    //double radius = 0.6;

    /* setup arguments:
    double radius, int max_length, int max_nb_seeds, 
    int avoid_u_turns, double max_angle, 
    double singularity_mask_radius,
    unsigned int random_seed, int parallel, int num_threads */
    engine.setup(radius, 1024, 32, 1, 90, 0, 0.1, 129873, 0, -1);
    //int rgn[] = {50, 100, 130};
    //engine.define_seed_region(3, rgn);
    engine.run();
    
    for (int i=0;i<engine.collection.size();i++) {
        engine.collection[i]->print();
        printf("\n");
    }
}
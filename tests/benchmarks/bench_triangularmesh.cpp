#include "engine.cpp"
#include "tictoc.h"
#include "readmesh.cpp"

int main()
{
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off("../data/test_finemesh2d.off", nv, nt, ver, tri);
    read_orient_txt("../data/test_finemesh2d.orient", nt, orient);
    TriangularMesh S;
    S.initialize(nv, nt, ver, tri, orient, 0, 1);
    double radius = 0.02;

    tic();
    S.identify_incident_triangles();
    toc("indicent triangles");
    S.identify_adjacent_triangles();
    toc("adjacent triangles");
    S.identify_triangle_close_neighborhood();
    toc("close neighborhood");
    S.create_faces();
    toc("faces");
    S.identify_boundary();
    toc("boundary");
    
    std::vector<int> focus, uturn, node;
    S.find_singularities(focus, uturn, node);
    toc("singularities");

    S.calculate_bounding_spheres();
    toc("bounding spheres");
    S.identify_triangle_extended_neighborhood(radius);
    toc("extended neighborhood");
    S.identify_triangle_extended_neighborhood_parallel(radius);
    toc("extended neighborhood parallel");
    
    printf("radius / min altitude = %.4f\n", radius/S.triangle_min_altitude);
    printf("%d (focus) + %d (u-turn) + %d (node) singularities\n", focus.size(), uturn.size(), node.size());
}
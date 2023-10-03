#include "engine.cpp"
#include "readmesh.cpp"

int main()
{
    char fname[] = "../data/test_mesh2d.off";
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off(fname, nv, nt, ver, tri);
    dummy_orient(nt, orient);
    TriangularMesh S;
    S.initialize(nv, nt, ver, tri, orient, 0, 0);
    double radius = 2.0;
    S.compute_neighborhood(radius, 0);

    // check t_neigh_tri
    FILE* file = fopen("test_mesh2d.neightri", "w");
    fprintf(file, "%.5f\n", radius);
    for (int i=0;i<S.nt;i++) {
        for (int j=0;j<S.t_nb_extneigh_tri[i];j++)
            fprintf(file, "%d ", S.t_extneigh_tri[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
}
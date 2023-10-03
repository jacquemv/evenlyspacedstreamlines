#include "engine.cpp"
#include "readmesh.cpp"

int main()
{
    char fname[] = "../data/test_mesh3d.off";
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off(fname, nv, nt, ver, tri);
    dummy_orient(nt, orient);

    TriangularMesh S;
    S.initialize(nv, nt, ver, tri, orient, 0, 0);

    FILE* file = fopen("test_mesh3d.nbound", "w");
    fprintf(file, "%.9f", 2 - S.euler_characteristic);
    fclose(file);

    file = fopen("test_mesh3d.radius", "w");
    fprintf(file, "%.9f", S.triangle_min_altitude);
    fclose(file);

    // check v_incident_tri
    file = fopen("test_mesh3d.inctri", "w");
    for (int i=0;i<S.nv;i++) {
        for (int j=0;j<S.v_nb_incident_tri[i];j++)
            fprintf(file, "%d ", S.v_incident_tri[S.v_maxtri*i+j]);
        fprintf(file, "\n");
    }
    fclose(file);

    // check t_adjacent_tri
    file = fopen("test_mesh3d.adjtri", "w");
    for (int i=0;i<S.nt;i++) {
        fprintf(file, "%d %d %d\n", S.t_adjacent_tri[3*i], 
                                    S.t_adjacent_tri[3*i+1], 
                                    S.t_adjacent_tri[3*i+2]);
    }
    fclose(file);

    // check is_boundary_vertex
    file = fopen("test_mesh3d.bound", "w");
    for (int i=0;i<S.nv;i++) if (S.is_boundary_vertex[i]) fprintf(file, "%d ", i);
    fclose(file);

    // check t_neigh_tri
    file = fopen("test_mesh3d.neightri", "w");
    for (int i=0;i<S.nt;i++) {
        for (int j=0;j<S.t_nb_neigh_tri[i];j++)
            fprintf(file, "%d ", S.t_neigh_tri[S.t_maxtri*i+j]);
        fprintf(file, "\n");
    }
    fclose(file);

    // check bounding sphere
    /*file = fopen("test_mesh3d.spheres", "w");
    for (int i=0;i<S.nt;i++) {
        double *c = S.bounding_sphere + 4*i;
        fprintf(file, "%.4f %.4f %.4f  %.4f\n", c[0], c[1], c[2], c[3]);
    }
    fclose(file);*/
}
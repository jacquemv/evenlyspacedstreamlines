#include "engine.cpp"
#include "readmesh.cpp"

const double EPS = 1e-10;

int main()
{
    int nv, nt;
    double *ver, *orient;
    int *tri;
    read_off("../data/test_mesh2d.off", nv, nt, ver, tri);
    read_orient_txt("../data/test_mesh2d.orient", nt, orient);
    TriangularMesh S;
    S.initialize(nv, nt, ver, tri, orient, 0, 0);

    // test that the transformation of coordinate between adjacent triangles are inverse of each other
    for (int i=0;i<S.nt;i++) {
        double ai[3], bi[3], cosa[3];
        bool cont[3];
        int ni[3];
        S.get_triangle_connections(i, ni, ai, bi, cont, cosa);
        for (int k=0;k<3;k++) {
            int j = ni[k];
            if (j == -1) continue;
            int nj[3];
            double aj[3], bj[3];
            S.get_triangle_connections(j, nj, aj, bj, cont, cosa);
            int m = (i == nj[1]) + 2*(i == nj[2]);
            assert( fabs(ai[k]*aj[m] - 1) < EPS && fabs(ai[k]*bj[m]+bi[k]) < EPS );
        }
    }
}
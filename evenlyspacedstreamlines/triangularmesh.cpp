#include <math.h>
#include "algebra.h"
#include "geometry.h"
#include "triangle.h"
#include "triangularmesh.h"

//-----------------------------------------------------------------------------
TriangularMesh::TriangularMesh()
{
    nv = nt = 0;
    ver = (double*) NULL;
    tri = (int*) NULL;
    orient = (double*) NULL;
    orthogonal = 0;

    face = (Triangle*) NULL;
    v_incident_tri = (int*) NULL;
    v_nb_incident_tri = (int*) NULL;
    t_adjacent_tri = (int*) NULL;
    t_nb_adjacent_tri = (int*) NULL;
    t_neigh_tri = (int*) NULL;
    t_nb_neigh_tri = (int*) NULL;
    is_boundary_vertex = (bool*) NULL;

    extneigh_allocated = 0;
    t_extneigh_tri = (int**) NULL; 
    t_nb_extneigh_tri = (int*) NULL;
    bounding_sphere = (double*) NULL;
}

//-----------------------------------------------------------------------------
int TriangularMesh::initialize(int nv_, int nt_, double* ver_, int* tri_, 
                                double* orient_, bool orthogonal_, 
                                bool allow_tweaking_orientation_)
{
    nv = nv_;
    nt = nt_;
    ver = ver_;
    tri = tri_;
    orient = orient_;
    orthogonal = orthogonal_;
    allow_tweaking_orientation = allow_tweaking_orientation_;

    t_extneigh_tri = (int**) NULL; 
    t_nb_extneigh_tri = (int*) NULL;
    bounding_sphere = (double*) NULL;
    extneigh_allocated = 0;

    return preprocess();
}

//-----------------------------------------------------------------------------
TriangularMesh::~TriangularMesh()
{
    delete [] face;
    delete [] v_incident_tri;
    delete [] v_nb_incident_tri;
    delete [] t_adjacent_tri;
    delete [] t_nb_adjacent_tri;
    delete [] t_neigh_tri;
    delete [] t_nb_neigh_tri;
    delete [] is_boundary_vertex;
    delete [] bounding_sphere;
    if (extneigh_allocated) {
        delete [] t_nb_extneigh_tri;
        for (int i=0;i<nt;i++) delete [] t_extneigh_tri[i];
    }
    delete [] t_extneigh_tri;
}

//-----------------------------------------------------------------------------
int TriangularMesh::get_triangle_connections(int idx, int* idx_adj, 
                                             double* a, double* b, 
                                             bool* is_continuous,
                                             double *cos_angle)
{
    idx_adj[0] = idx_adj[1] = idx_adj[2] = -1;
    a[0] = a[1] = a[2] = b[0] = b[1] = b[2] = 0;
    int n = t_nb_adjacent_tri[idx];
    for (int k=0;k<n;k++) {
        int i = t_adjacent_tri[3*idx+k];
        double s_crit = face[idx].s_crit;
        double s0 = face[i].vertex_coordinate(face[idx].idx_origin);
        double sc = face[i].vertex_coordinate(face[idx].idx_apex);
        double s1 = face[i].vertex_coordinate(face[idx].idx_base);
        int f = face[idx].common_edge(tri+3*i);
        idx_adj[f] = i;
        int i1=0, i2=0; // i1 -> i2 = common edge
        int i3=0; // i3 = other vertex in triangle idx
        int i4=0; // i4 = other vertex in triangle i

        // s' = a * s + b
        if (f == 0) { // commoen edge = origin-base
            // s = 0 -> s' = s0; s = 1 -> s' = s1
            b[f] = s0;
            a[f] = s1 - s0;
            i1 = face[idx].idx_origin;
            i2 = face[idx].idx_base;
            i3 = face[idx].idx_apex;
        }
        if (f == 1) { // commoen edge = origin-apex
            // s = 0 -> s' = s0; s = s_crit -> s' = sc
            b[f] = s0;
            a[f] = (sc - s0)/s_crit;
            i1 = face[idx].idx_origin;
            i2 = face[idx].idx_apex;
            i3 = face[idx].idx_base;
        }
        if (f == 2) { // commoen edge = base-apex
            // s = 1 -> s' = s1; s = s_crit -> s' = sc
            a[f] = (s1 - sc)/(1 - s_crit);
            b[f] = s1 - a[f];
            i1 = face[idx].idx_base;
            i2 = face[idx].idx_apex;
            i3 = face[idx].idx_origin;
        }

        i4  = tri[3*i] + tri[3*i+1] + tri[3*i+2] - i1 - i2; // remaining index
        double common_edge[3];
        vdiff(common_edge, ver+3*i1, ver+3*i2);

        double outward3[3], outward4[3]; // away from the common edge
        vdiff(outward3, ver+3*i1, ver+3*i3);
        vrmcomp(outward3, common_edge);
        bool b3 = vdot(face[idx].orient, outward3) < 0;

        vdiff(outward4, ver+3*i1, ver+3*i4);
        vrmcomp(outward4, common_edge);
        bool b4 = vdot(face[i].orient, outward4) < 0;

        cos_angle[f] = vdot(face[i].orient, face[idx].orient) 
                          / face[i].height / face[idx].height;
        if (b3 == b4)
            cos_angle[f] = -cos_angle[f];
        
        is_continuous[f] = b3 ^ b4 ^ face[i].reverse_orient 
                                   ^ face[idx].reverse_orient;
    }
    return n;
}

//-----------------------------------------------------------------------------
void TriangularMesh::get_triangle_parameters(int idx, double& s_crit, 
                                         double& height)
{
    s_crit = face[idx].s_crit;
    height = face[idx].height;
}

//-----------------------------------------------------------------------------
int TriangularMesh::preprocess()
{
    identify_incident_triangles();
    if (identify_adjacent_triangles() < 0) {
        return -1;
    }
    identify_boundary();
    identify_triangle_close_neighborhood();
    create_faces();

    triangle_min_altitude = calculate_triangle_min_altitude();
    triangle_max_base = calculate_triangle_max_base();
    euler_characteristic = calculate_euler_characteristic();
    neighborhood_mean_size = 0.0;
    return 0;
}

//-----------------------------------------------------------------------------
void TriangularMesh::create_faces()
{
    face = new Triangle [nt];
    for (int i=0;i<nt;i++)
        face[i] = Triangle(ver, tri+3*i, orient+3*i, orthogonal, 
                           allow_tweaking_orientation);
}

//-----------------------------------------------------------------------------
void TriangularMesh::identify_incident_triangles()
{
    v_nb_incident_tri = new int [nv];
    for (int i=0;i<nv;i++) v_nb_incident_tri[i] = 0;
    for (int i=0;i<3*nt;i++) v_nb_incident_tri[tri[i]]++;

    v_maxtri = 3;
    for (int i=0;i<nv;i++)
        if (v_nb_incident_tri[i] > v_maxtri)
            v_maxtri = v_nb_incident_tri[i];

    v_incident_tri = new int [nv*v_maxtri];
    for (int i=0;i<nv;i++) v_nb_incident_tri[i] = 0;
    for (int i=0;i<3*nt;i++) {
        v_incident_tri[v_maxtri*tri[i] + (v_nb_incident_tri[tri[i]]++)] = i/3;
    }
}

//-----------------------------------------------------------------------------
void TriangularMesh::find_singularities(std::vector<int> &focus, 
                                        std::vector<int> &uturn,
                                        std::vector<int> &node)
{
    focus.clear(); uturn.clear(); node.clear();
    for (int i=0;i<nv;i++) {
        // # ways out of the vertex following an orienttion vector
        int ways_out = 0; 
        int nb_tri = v_nb_incident_tri[i];
        for (int j=0;j<nb_tri;j++) {
            int idtri = v_incident_tri[v_maxtri*i+j];
            
            assert( i == tri[3*idtri] || i == tri[3*idtri+1] 
                 || i == tri[3*idtri+2] );
            assert( face[idtri].idx_apex == tri[3*idtri] 
                 || face[idtri].idx_apex == tri[3*idtri+1] 
                 || face[idtri].idx_apex == tri[3*idtri+2] );

            ways_out += (i == face[idtri].idx_apex);
        }
        if (ways_out == 0) {
            focus.push_back(i);
        } else if (ways_out == 1) {
            if (!is_boundary_vertex[i])
                uturn.push_back(i);
        } else if (ways_out == nb_tri) {
            if (!is_boundary_vertex[i]) 
                node.push_back(i);
        }
    }
}

//-----------------------------------------------------------------------------
inline int count_identical(int* a, int* b)
{
    return (a[0] == b[0]) + (a[0] == b[1]) + (a[0] == b[2])
        +  (a[1] == b[0]) + (a[1] == b[1]) + (a[1] == b[2])
        +  (a[2] == b[0]) + (a[2] == b[1]) + (a[2] == b[2]);
}

//-----------------------------------------------------------------------------
inline int count_matches(int i, int*a, int* b)
{
    return (a[0] == i) + (a[1] == i) + (a[2] == i) 
         + (b[0] == i) + (b[1] == i) + (b[2] == i);
}

//-----------------------------------------------------------------------------
int TriangularMesh::identify_adjacent_triangles()
{
    t_nb_adjacent_tri = new int [nt];
    for (int i=0;i<nt;i++)
        t_nb_adjacent_tri[i] = 0;
    t_adjacent_tri = new int [3*nt];
    for (int i=0;i<3*nt;i++)
        t_adjacent_tri[i] = -1;

    for (int i=0;i<3*nt;i++) {
        int i1 = i / 3;
        int k = tri[i];
        for (int j=0;j<v_nb_incident_tri[k];j++) {
            int i2 = v_incident_tri[v_maxtri*k+j];
            if (count_identical(tri+3*i1, tri+3*i2) == 2) {
                if (   (t_adjacent_tri[3*i1] == i2) 
                    || (t_adjacent_tri[3*i1+1] == i2) 
                    || (t_adjacent_tri[3*i1+2] == i2))
                    continue;
                // otherwise the sirface is not a manifold (could be an error)
                if (t_nb_adjacent_tri[i1] < 3)
                    t_adjacent_tri[3*i1+(t_nb_adjacent_tri[i1]++)] = i2;
                else
                    return -1;
            }
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------
void TriangularMesh::identify_boundary()
{
    is_boundary_vertex = new bool [nv];
    for (int i=0;i<nv;i++)
        is_boundary_vertex[i] = 0;
    for (int i=0;i<nt;i++) {
        if (t_nb_adjacent_tri[i] == 3)
            continue;
        if (t_nb_adjacent_tri[i] <= 1) {
            for (int k=0;k<3;k++) 
                is_boundary_vertex[tri[3*i+k]] = 1;
        } else {
            int i0 = t_adjacent_tri[3*i], i1 = t_adjacent_tri[3*i+1];
            for (int k=0;k<3;k++) {
                int idx = tri[3*i+k];
                if (count_matches(idx, tri+3*i0, tri+3*i1) == 1)
                    is_boundary_vertex[idx] = 1;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void TriangularMesh::identify_triangle_close_neighborhood()
{
    t_maxtri = 3*v_maxtri - 3;
    t_nb_neigh_tri = new int [nt];
    t_neigh_tri = new int [nt*t_maxtri];
    int t_maxtri_opt = 0;
    for (int i=0;i<nt;i++) {
        // start with adjacent faces
        t_nb_neigh_tri[i] = t_nb_adjacent_tri[i];
        t_neigh_tri[i*t_maxtri] = t_adjacent_tri[3*i];
        t_neigh_tri[i*t_maxtri+1] = t_adjacent_tri[3*i+1];
        t_neigh_tri[i*t_maxtri+2] = t_adjacent_tri[3*i+2];
        // add other faces, but not the current face
        for (int k=0;k<3;k++) {
            int v = tri[3*i+k];
            for (int n=0;n<v_nb_incident_tri[v];n++) {
                int t = v_incident_tri[v*v_maxtri+n];
                if (   (t == i) 
                    || (t == t_adjacent_tri[3*i]) 
                    || (t == t_adjacent_tri[3*i+1]) 
                    || (t == t_adjacent_tri[3*i+2]))
                    continue;
                t_neigh_tri[i*t_maxtri+(t_nb_neigh_tri[i]++)] = t;
            }
        }
        if (t_maxtri_opt < t_nb_neigh_tri[i])
            t_maxtri_opt = t_nb_neigh_tri[i];
    }
    // compress table
    for (int i=0;i<nt;i++)
        for (int j=0;j<t_nb_neigh_tri[i];j++)
            t_neigh_tri[i*t_maxtri_opt+j] = t_neigh_tri[i*t_maxtri+j];
    t_maxtri = t_maxtri_opt;
}

//-----------------------------------------------------------------------------
void TriangularMesh::calculate_bounding_spheres()
{
    bounding_sphere = new double [4*nt];
    for (int i=0;i<nt;i++) {
        double *c = bounding_sphere + 4*i;
        face[i].bounding_sphere(c, c[3]);
    }
}

//-----------------------------------------------------------------------------
bool TriangularMesh::triangles_close_enough(int i, int j, double radius)
{
    double radius2 = radius*radius;
    // check bounding spheres
    double *Gi = bounding_sphere + 4*i;
    double *Gj = bounding_sphere + 4*j;
    double Ri = Gi[3], Rj = Gj[3];
    double R = vdist(Gi, Gj);
    if (R > Ri + Rj + radius)
        return 0;
    if (R < radius) return 1;

    double *a1, *b1, *c1, *a2, *b2, *c2;
    a1 = ver + 3*tri[3*i];
    b1 = ver + 3*tri[3*i+1];
    c1 = ver + 3*tri[3*i+2];
    a2 = ver + 3*tri[3*j];
    b2 = ver + 3*tri[3*j+1];
    c2 = ver + 3*tri[3*j+2];

    if (segment_distance2(a1, b1, a2, b2) < radius2) return 1;
    if (segment_distance2(a1, b1, c2, b2) < radius2) return 1;
    if (segment_distance2(a1, b1, c2, a2) < radius2) return 1;
    if (segment_distance2(a1, c1, a2, b2) < radius2) return 1;
    if (segment_distance2(a1, c1, c2, b2) < radius2) return 1;
    if (segment_distance2(a1, c1, c2, a2) < radius2) return 1;
    if (segment_distance2(b1, c1, a2, b2) < radius2) return 1;
    if (segment_distance2(b1, c1, c2, b2) < radius2) return 1;
    if (segment_distance2(b1, c1, c2, a2) < radius2) return 1;
    return 0;
}

//-----------------------------------------------------------------------------
void TriangularMesh::identify_triangle_extended_neighborhood(double radius)
{
    t_extneigh_tri = new int* [nt];
    t_nb_extneigh_tri = new int [nt];
    extneigh_allocated = 1;

    int* tag = new int [nt];
    for (int i=0;i<nt;i++) tag[i] = -1;

    int prev_nb, new_nb, total_nb;
    int* prev_list = new int [nt];
    int* new_list = new int [nt];
    int* total_list = new int [nt];

    int count_true = 0, count_false = 0;

    for (int i=0;i<nt;i++) {
        prev_nb = total_nb = t_nb_neigh_tri[i];
        tag[i] = i;
        for (int j=0;j<prev_nb;j++) {
            int idx = t_neigh_tri[i*t_maxtri+j];
            prev_list[j] = total_list[j] = idx;
            tag[idx] = i;
        }

        while (1) {
            // find list of next order neighbors
            new_nb = 0;
            for (int j=0;j<prev_nb;j++) {
                int idx = prev_list[j];
                for (int k=0;k<t_nb_neigh_tri[idx];k++) {
                    int idx_neigh = t_neigh_tri[idx*t_maxtri+k];
                    if (tag[idx_neigh] != i) {
                        tag[idx_neigh] = i;
                        new_list[new_nb++] = idx_neigh;
                    }
                }
            }

            int cnt = 0;
            for (int j=0;j<new_nb;j++) {
                if (triangles_close_enough(i, new_list[j], radius)) {
                    total_list[total_nb++] = new_list[j];
                    cnt++;
                    count_true++;
                } else count_false++;
            }
            if (cnt == 0) break;

            // swap lists
            int *tmp, tmpi;
            tmp = new_list;
            new_list = prev_list;
            prev_list = tmp;
            tmpi = new_nb;
            new_nb = prev_nb;
            prev_nb = tmpi;

        }

        t_extneigh_tri[i] = new int [total_nb];
        t_nb_extneigh_tri[i] = total_nb;
        for (int j=0;j<total_nb;j++) t_extneigh_tri[i][j] = total_list[j];
    }

    delete [] tag;
    delete [] prev_list;
    delete [] new_list;
    delete [] total_list;
}

//-----------------------------------------------------------------------------
void TriangularMesh::identify_triangle_extended_neighborhood_parallel(double radius)
{
    t_extneigh_tri = new int* [nt];
    t_nb_extneigh_tri = new int [nt];
    extneigh_allocated = 1;
    int count_true = 0, count_false = 0;

    #pragma omp parallel
    {
        int* tag = new int [nt];
        for (int i=0;i<nt;i++) tag[i] = -1;

        int prev_nb, new_nb, total_nb;
        int* prev_list = new int [nt];
        int* new_list = new int [nt];
        int* total_list = new int [nt];

        #pragma omp for schedule(static,64) reduction(+:count_true) reduction(+:count_false)
        for (int i=0;i<nt;i++) {
            prev_nb = total_nb = t_nb_neigh_tri[i];
            tag[i] = i;
            for (int j=0;j<prev_nb;j++) {
                int idx = t_neigh_tri[i*t_maxtri+j];
                prev_list[j] = total_list[j] = idx;
                tag[idx] = i;
            }

            while (1) {
                // find list of next order neighbors
                new_nb = 0;
                for (int j=0;j<prev_nb;j++) {
                    int idx = prev_list[j];
                    for (int k=0;k<t_nb_neigh_tri[idx];k++) {
                        int idx_neigh = t_neigh_tri[idx*t_maxtri+k];
                        if (tag[idx_neigh] != i) {
                            tag[idx_neigh] = i;
                            new_list[new_nb++] = idx_neigh;
                        }
                    }
                }

                int cnt = 0;
                for (int j=0;j<new_nb;j++) {
                    if (triangles_close_enough(i, new_list[j], radius)) {
                        total_list[total_nb++] = new_list[j];
                        cnt++;
                        count_true++;
                    } else count_false++;
                }
                if (cnt == 0) break;

                // swap lists
                int *tmp, tmpi;
                tmp = new_list;
                new_list = prev_list;
                prev_list = tmp;
                tmpi = new_nb;
                new_nb = prev_nb;
                prev_nb = tmpi;

            }

            t_extneigh_tri[i] = new int [total_nb];
            t_nb_extneigh_tri[i] = total_nb;
            for (int j=0;j<total_nb;j++) t_extneigh_tri[i][j] = total_list[j];
        }
        delete [] tag;
        delete [] prev_list;
        delete [] new_list;
        delete [] total_list;
    } // end parallel region
}

//-----------------------------------------------------------------------------
void TriangularMesh::default_extended_neighborhood()
{
    t_extneigh_tri = new int* [nt];
    t_nb_extneigh_tri = t_nb_neigh_tri;
    extneigh_allocated = 0;
    for (int i=0;i<nt;i++) {
        t_extneigh_tri[i] = t_neigh_tri + t_maxtri * i;
    }
}

//-----------------------------------------------------------------------------
void TriangularMesh::compute_neighborhood(double radius, int parallel)
{

    if (radius <= triangle_min_altitude) {
        default_extended_neighborhood();
    } else {
        calculate_bounding_spheres();
        if (parallel)
            identify_triangle_extended_neighborhood_parallel(radius);
        else
            identify_triangle_extended_neighborhood(radius);
    }
    neighborhood_mean_size = calculate_neighborhood_mean_size();
}

//-----------------------------------------------------------------------------
inline void TriangularMesh::get_neighborhood(int idx, int* &idx_list, 
                                             int& list_len)
{
    list_len = t_nb_extneigh_tri[idx];
    idx_list = t_extneigh_tri[idx];
}

//-----------------------------------------------------------------------------
 inline void TriangularMesh::get_incident_triangles(int idver, int* &idx_list, 
                                                    int& list_len)
 {
    list_len = v_nb_incident_tri[idver];
    idx_list = v_incident_tri + v_maxtri*idver;
 }

//-----------------------------------------------------------------------------
 double TriangularMesh::calculate_neighborhood_mean_size()
 {
    int sum = 0;
    if (t_extneigh_tri == NULL)
        for (int i=0;i<nt;i++)
            sum += t_nb_neigh_tri[i];
    else
        for (int i=0;i<nt;i++)
            sum += t_nb_extneigh_tri[i];
    return ((double) sum) / nt;
 }

//-----------------------------------------------------------------------------
inline void TriangularMesh::compute_mask(int idx, double* P, double* Q, 
                                         double radius, 
                                         double& s1, double& s2)
{
    face[idx].coordinate_mask_cylinder(P, Q, radius, s1, s2);
}

//-----------------------------------------------------------------------------
inline double TriangularMesh::shrink_factor(int idx)
{
    return face[idx].shrink;
}

//-----------------------------------------------------------------------------
inline void TriangularMesh::get_segment(int idx, double s, 
                                        double* P, double* Q)
{
    face[idx].segment(s, P, Q);
}

//-----------------------------------------------------------------------------
inline void TriangularMesh::get_partial_segment(int idx, double s, 
                                    double t1, double t2, double* P, double* Q)
{
    double A[3], B[3];
    face[idx].segment(s, A, B);
    vlincomb(P, 1-t1, A, t1, B);
    vlincomb(Q, 1-t2, A, t2, B);
}

//-----------------------------------------------------------------------------
double TriangularMesh::calculate_triangle_min_altitude()
{
    // compute triangle altitudes as 2*area/base
    double R2 = 9e99;
    for (int i=0;i<nt;i++) {
        double *A = ver+3*tri[3*i], 
               *B = ver+3*tri[3*i+1], 
               *C = ver+3*tri[3*i+2];
        double v[3], w[3], n[3];
        vdiff(v, A, B);
        vdiff(w, A, C);
        vcross(n, v, w);
        double h2, A2 = vnorm2(n);
        if (A2 <= 0) return 0;
        h2 = A2/vnorm2(v); if (h2 < R2) R2 = h2;
        h2 = A2/vnorm2(w); if (h2 < R2) R2 = h2;
        vdiff(v, C, B);
        h2 = A2/vnorm2(v); if (h2 < R2) R2 = h2;
    }
    return sqrt(R2);
}

//-----------------------------------------------------------------------------
double TriangularMesh::calculate_triangle_max_base()
{
    double max_base2 = 0;
    for (int i=0;i<nt;i++)
        if (max_base2 < face[i].b2) max_base2 = face[i].b2;
    return sqrt(max_base2);
}

//-----------------------------------------------------------------------------
int TriangularMesh::calculate_euler_characteristic()
{
    int count[] = {0, 0, 0, 0};
    for (int i=0;i<nt;i++)
        count[t_nb_adjacent_tri[i]]++;
    
    int ne = (6*count[0]+5*count[1]+4*count[2]+3*count[3])/2;
    int chi = nv - ne + nt;
    return chi;
}

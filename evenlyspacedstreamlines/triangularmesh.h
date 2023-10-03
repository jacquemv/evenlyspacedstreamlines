#ifndef TRIANGULARMESH_H_
#define TRIANGULARMESH_H_

#include <vector>

//-----------------------------------------------------------------------------
class TriangularMesh {
public:
    // triangular mesh
    int nv, nt; // number of vertices and triangles
    double* ver; // 3*nv vertex positions
    int* tri; // 3*nt vertex indices of triangles
    double* orient; // 3*nt orientation vectors
    bool orthogonal;// rotate orientation by 90 degrees
    bool allow_tweaking_orientation; // add noise to the orientation  
                                     // if parallel to an edge

    TriangularMesh();
    ~TriangularMesh();
    
    // returns 0 if no error occured
    int initialize(int nv_, int nt_, double* ver_, int* tri_, 
                   double* orient_, bool orthogonal_, 
                   bool allow_tweaking_orientation_);

    // mesh global measurements
    double triangle_min_altitude;
    double triangle_max_base;
    int euler_characteristic;
    double neighborhood_mean_size;
    
    // get information about triangle idx
    void get_triangle_parameters(int idx, double& s_crit, double& height);
    // find the mapping between local coordinates of adjacent triangles
    // as well as the cosine of the angle between orientation vectors
    int get_triangle_connections(int idx, int* idx_adj, double* a, double* b, 
                                 bool* is_continuous, double *cos_angle);
    // scaling factor from xyz coordinates to s coordinate
    inline double shrink_factor(int idx);

    // find the list of vertex indices that are singularities of the orientation
    // field. focus = no way out of the vertex; uturn = one single way out;
    // node = a way out in every incident triangle (attractor/repellor)
    void find_singularities(std::vector<int> &focus, std::vector<int> &uturn,
                            std::vector<int> &node);

    // list of triangles at a distance <= radius from a triangle
    void compute_neighborhood(double radius, int parallel);
    inline void get_neighborhood(int idx, int* &idx_list, int& list_len);
    // list of triangles incident to a vertex
    inline void get_incident_triangles(int idver, int* &idx_list, int& list_len);

    // xyz coordinates of a segment at coordinate s on triangle idx
    inline void get_segment(int idx, double s, double* P, double* Q);
    inline void get_partial_segment(int idx, double s, 
                                    double t1, double t2, double* P, double* Q);
    
    // mask of s coordinates cast by the segment PQ at distance radius
    inline void compute_mask(int idx, double* P, double* Q, double radius, 
                             double& s1, double& s2);

//private:
    // triangles incident to a vertex
    int* v_incident_tri;
    int* v_nb_incident_tri;
    int v_maxtri;

    // triangles adjacent to a triangle
    int* t_adjacent_tri;
    int* t_nb_adjacent_tri;

    // triangles in the neighborhood of a triangle
    int* t_neigh_tri; 
    int* t_nb_neigh_tri;
    int t_maxtri;

    // triangles in the extended neighborhood of a triangle
    int** t_extneigh_tri; 
    int* t_nb_extneigh_tri;
    double* bounding_sphere;
    bool extneigh_allocated;

    // boundary vertices
    bool* is_boundary_vertex;

    // triangle objects
    Triangle* face;

    // preprocessing
    int preprocess(); // execute all the functions below
    void identify_incident_triangles();
    int identify_adjacent_triangles();
    void identify_boundary();
    void identify_triangle_close_neighborhood();
    void create_faces();

    double calculate_triangle_min_altitude();
    double calculate_triangle_max_base();
    int calculate_euler_characteristic();
    double calculate_neighborhood_mean_size();

    // extended neighborhood
    void default_extended_neighborhood();
    void identify_triangle_extended_neighborhood(double radius);
    void identify_triangle_extended_neighborhood_parallel(double radius);
    void calculate_bounding_spheres();
    bool triangles_close_enough(int i, int j, double radius);
};

#endif
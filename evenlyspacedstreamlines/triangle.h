#ifndef TRIANGLE_H_
#define TRIANGLE_H_

//-----------------------------------------------------------------------------
class Triangle {
public:
    double* ver; // vertex positions
    int idx_origin, idx_base, idx_apex; // vertex indices

    // Local basis on the triangle:
    // the indices idx_origin, idx_base, idx_apex have the same orientation
    // as the original triangle

    // In the basis (base_vector, orient), a point x within the triangle is
    // expressed as
    //   x = s * base_vector + t * orient
    // such that
    //   0 < s < 1
    //   0 < t < s/s_crit             if 0 < s <= s_crit
    //   0 < t < (1-s)/(1-s_crit)     if s_crit <= s <`1
    // in particular, the vertices (idx_origin, idx_base, idx_apex) are
    // mapped to the coordinates (0, 0), (1, 0), (s_crit, 1)
    // and in that system, orient becomes vertical

    double origin[3]; // position of idx_origin
    double base_vector[3]; // vector from idx_origin to idx_base
    double orient[3]; // orientation vector
    double s_crit, height;
    double shrink; // shrink factor in the transverse direction from xyz to s

    Triangle() {}
    Triangle(double* ver_, int* idx, double* orient, bool flip_orient, 
             bool allow_tweaking_orientation);

    // map from the local basis to the 3D space
    void local_to_global(double s, double t, double* P);

    // project from the 3D space to local basis
    void global_to_local(double* P, double& s, double& t);

    // returns the s coordinate of one of the three vertex of the triangle
    double vertex_coordinate(int idx);

    // returns 1 if (s, t) is in the triangle and 0 otherwise
    inline bool contains(double s, double t);

    // determines which edge the triangle has in common with another triangle
    // returns 0 (origin-base), 1 (origin-apex), 2 (base-apex), -1 (none)
    int common_edge(int* idx);

    // segment PQ on the triangle parallel to the orientation vector
    // at coordinate s along base_vector
    void segment(double s, double* P, double* Q);

    // length of that segment
    double segment_length(double s);

    // smallest enclosing sphere
    void bounding_sphere(double *center, double &radius);

    // range (s_min, s_max) of s coordinates in the triangle that are at a 
    // distance less than R from a segment PQ
    int coordinate_mask_cylinder(double* P, double* Q, double R, 
                                 double& s_min, double& s_max);

//private:
    double b2, h2, bh; // scalar products
    bool reverse_orient;

    // intersection between a segment in the local coordinates and a sphere
    int intersection_segment_sphere(double* P, double R, // center and radius
                                    double s1, double t1, // point 1
                                    double s2, double t2, // point 2
                                    double& si1, double& ti1, // intersection 1
                                    double& si2, double& ti2 // intersection 2
                                   );

    int coordinate_mask_quadratics(EllipticQuadraticForm qf,
                                   double& s_min, double& s_max);
    
    int coordinate_mask_quadratics_constrained(EllipticQuadraticForm qf,
                                double* line, double& s_min, double& s_max);
    
    void linear_constraint(double* P, double* Q, double* line);
    
    int coordinate_mask_sphere(double* P, double R, 
                               double& s_min, double& s_max);
                            
    int coordinate_mask_infinite_cylinder(double* P, double* Q, double R, 
                                 double& s_min, double& s_max);
};

#endif
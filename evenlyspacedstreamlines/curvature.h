#ifndef CURVATURE_H_
#define CURVATURE_H_

class Curvature {
public:
    Curvature();
    Curvature(int max_nb_segments_);
    ~Curvature();

    void initialize(int max_nb_segments_);
    void reset();
    void append_segment(double* p, double* q);
    void low_curvature_interval(double radius, int &start, int &len);

//private:
    int nb_segments, max_nb_segments;
    double* points; // size 6*nb_segments
    int *constraint_lo, *constraint_hi;
    bool calculate_constraint(double radius);
};

#endif
#ifndef ELLIPTICQUADRATICFORM_H_
#define ELLIPTICQUADRATICFORM_H_

//-----------------------------------------------------------------------------
class EllipticQuadraticForm {
public:
    /*  The quadratic form
            Q(x, y) = a x^2 + b y^2 + 2 c x y + 2 d x + 2 e y + f
        defines a region Q(x, y) <= 0.

        We assume that Q(x, y) has a lower bound (as in distance problems)
        which leads to the following cases:
        - Ellipse: a > 0, b > 0, a b-c^2 > 0
        - Infinite strip: a >= 0, b >= 0, a b = c^2, b d = c e
        - Empty set
    */
    double a, b, c, d, e, f; 
    double delta, determinant;

    // create quadratic form
    EllipticQuadraticForm() {}
    EllipticQuadraticForm(double a_, double b_, double c_, 
                          double d_, double e_, double f_);

    // create an ellipse with center (cx, cy), two axis long_axis and 
    // short_axis, rotated with an angle in radians
    void ellipse(double xc, double yc, double long_axis, double short_axis, 
                 double angle);
    
    // create an infinite strip with a width ± halfwidth on both sides of a line
    // passing through the center (xc, yc) with an angle in radians
    void infinite_strip(double xc, double yc, double halfwidth, double angle);

    // compute Q(x, y)
    double evaluate(double x, double y);

    // compute the global minimum of the quadratic form
    double minimum();

    // find the center of the ellipse; for an infinite strip pick one point
    // on the line of minima
    // returns 1 if successful, and 0 otherwise
    int center(double& xc, double& yc);

    // find the two axes of the ellipse and returns 1 if successful,
    // and 0 otherwise
    int semiaxis(double& long_axis, double& short_axis);

    // find the half-width of an infinite strip or the short axis of the ellipse 
    // returns 1 if successful, and 0 if the set is empty
    int halfwidth(double& hwidth);

    // find the angle of the ellipse or the infinite strip (between -pi/2 and 
    // +pi/2) and returns 1 if successful and 0 otherwise
    int angle(double& theta);

    // find the extereme points (horizontal and vertical tangent points)
    // points (x1, y1) and (x2, y2) are ordered with respect to x (vertical 
    // tangents) or y (horizontal tangents)
    // returns 2 if successful, and 0 otherwise
    int extreme_points_xaxis(double& x1, double& y1, 
                             double& x2, double& y2);
    int extreme_points_yaxis(double& x1, double& y1, 
                             double& x2, double& y2);

    // compute the intersection with a line segment
    // returns the number of intersection points (0, 1, or 2)
    // tangent points are not included; in the case of
    int intersection_segment(double x1, double y1, 
                             double x2, double y2,
                             double& xi1, double& yi1,
                             double& xi2, double& yi2);
};

#endif
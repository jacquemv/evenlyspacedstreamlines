#include "quadraticforms.h"
#include <math.h>

//-----------------------------------------------------------------------------
EllipticQuadraticForm::EllipticQuadraticForm(double a_, double b_, double c_, 
                                             double d_, double e_, double f_)
{
    a = a_; b = b_; c = c_; d = d_; e = e_; f = f_;
    // delta > 0: ellipse ; delta = 0: infinite strip
    delta = a*b - c*c;
    if (delta < 0) delta = 0; // in case of rounding error
    // determinant = det([a c d; c b e; d e f]) <= 0 otherwise the set is empty
    determinant = f*delta + 2*c*d*e - a*e*e - b*d*d; 
}

//-----------------------------------------------------------------------------
void EllipticQuadraticForm::ellipse(double xc, double yc, 
                                    double long_axis, double short_axis, 
                                    double angle)
{
    double A2 = long_axis*long_axis;
    double B2 = short_axis*short_axis;
    double sina = sin(angle), cosa = cos(angle);
    double sin2a = sina*sina, cos2a = cosa*cosa;
    a = cos2a/A2 + sin2a/B2;
    b = sin2a/A2 + cos2a/B2;
    c = -sina * cosa * (1/B2 - 1/A2);
    d = -xc*a - yc*c;
    e = -yc*b - xc*c;
    f = -1 + a*xc*xc + b*yc*yc + 2*xc*yc*c;
    delta = a*b - c*c;
    determinant = f*delta + 2*c*d*e - a*e*e - b*d*d; 
}

//-----------------------------------------------------------------------------
void EllipticQuadraticForm::infinite_strip(double xc, double yc, 
                                          double halfwidth, double angle)
{
    double sina = sin(angle), cosa = cos(angle);
    double sin2a = sina*sina, cos2a = cosa*cosa;
    a = sin2a;
    b = cos2a;
    c = - sina * cosa;
    d = -xc*a - yc*c;
    e = -yc*b - xc*c;
    f = xc*xc*a + yc*yc*b + 2*xc*yc*c - halfwidth*halfwidth;
    delta = determinant = 0;
}

//-----------------------------------------------------------------------------
double EllipticQuadraticForm::evaluate(double x, double y)
{
    return a*x*x + b*y*y + 2*c*x*y + 2*d*x + 2*e*y + f;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::center(double& xc, double& yc)
{
    if (delta == 0) {
        if (a != 0) {
            xc = -d/a;
            yc = 0;
        } else {
            xc = 0;
            yc = -e/b;
        }
        return 1;
    }
    xc = (e*c-b*d)/delta;
    yc = (c*d-a*e)/delta;
    return 1;
}

//-----------------------------------------------------------------------------
double EllipticQuadraticForm::minimum()
{
    if (delta > 0)
        return determinant / delta; // true if a unique minimum exists
    if (a != 0) // the minimum is a line
        return f - d*d/a;
    return f - e*e/b;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::angle(double& theta)
{
    if (c != 0) {
        theta = atan((b-a - sqrt((a-b)*(a-b) + 4*c*c))/c/2);
        return 1;
    } else {
        theta = a < b ? 0 : 2*atan(1);
        return 1;
    }
    return 0;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::semiaxis(double& long_axis, double& short_axis)
{
    if ((determinant >=0) || (delta == 0)) return 0;
    double u = sqrt((a-b)*(a-b) + 4*c*c);
    long_axis = sqrt(-determinant/2 * (a+b + u))/delta;
    short_axis = sqrt(-determinant/2 * (a+b - u))/delta;
    return 1;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::halfwidth(double& hwidth)
{
    if (delta == 0) {
        if (a != 0) {
            double u = d*d - a*f;
            if (u <= 0) return 0;
            hwidth = sqrt(u/(a*a+c*c));
            return 1;
        } else {
            double u = e*e - b*f;
            if (u <= 0) return 0;
            hwidth = sqrt(u/(b*b+c*c));
            return 1;
        }
    }
    double dummy;
    return semiaxis(dummy, hwidth);
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::extreme_points_xaxis(double& x1, double& y1, 
                                                double& x2, double& y2)
{
    if (delta == 0) return 0;
    double u = d*b-c*e;
    double discr = u*u + delta*(e*e-b*f);
    if (discr <= 0) return 0;
    double sqrt_discr = sqrt(discr);
    x1 = (-u - sqrt_discr)/delta;
    x2 = (-u + sqrt_discr)/delta;
    y1 = -(c*x1+e)/b;
    y2 = -(c*x2+e)/b;
    return 2;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::extreme_points_yaxis(double& x1, double& y1, 
                                                double& x2, double& y2)
{
    if (delta == 0) return 0;
    double u = e*a-c*d;
    double discr = u*u + delta*(d*d-a*f);
    if (discr <= 0) return 0;
    double sqrt_discr = sqrt(discr);
    y1 = (-u - sqrt_discr)/delta;
    y2 = (-u + sqrt_discr)/delta;
    x1 = -(c*y1+d)/a;
    x2 = -(c*y2+d)/a;
    return 2;
}

//-----------------------------------------------------------------------------
int EllipticQuadraticForm::intersection_segment(double x1, double y1, 
                                                double x2, double y2,
                                                double& xi1, double& yi1, 
                                                double& xi2, double& yi2)
{
    // sort segment extremities
    if (x2 < x1) {
        double tmp;
        tmp = x1; x1 = x2; x2 = tmp;
        tmp = y1; y1 = y2; y2 = tmp;
    }

    // INFINITE STRIP

    if (delta == 0) {

        // horizontal strip
        if (a == 0) {
            if (y1 == y2) return 0;
            double w = e*e - b*f;
            if (w <= 0) return 0;
            w = sqrt(w);
            // x = p y + q
            int n = 0;
            double y = (-e-w)/b;
            double p = (x2-x1)/(y2-y1);
            double q = x1 - p*y1;
            double x = p*y + q;
            if ((x1 < x) && (x < x2)) {
                xi1 = x;
                yi1 = y;
                n = 1;
            }
            y += 2*w/b;
            x = x1 + p*(y-y1);
            if ((x1 < x) && (x < x2)) {
                if (n == 0) {
                    xi1 = x;
                    yi1 = y;
                } else {
                    xi2 = x;
                    yi2 = y;
                }
                n++;
            }
            return n;
        }

        double w = d*d - a*f;
        if (w <= 0) return 0;
        w = sqrt(w);

        // Vertical segment
        
        if (x1 == x2) {
            // x = x1
            if (c == 0) { // vertical strip
                return 0;
            }
            if (y2 < y1) {
                double tmp = y1; y1 = y2; y2 = tmp;
            }
            double y;
            int n = 0;
            y = (-a*x1 - d - w)/c;
            if ((y1 < y) && (y < y2)) {
                xi1 = x1;
                yi1 = y;
                n = 1;
            }
            y += 2*w/c;
            if ((y1 < y) && (y < y2)) {
                if (n == 0) {
                    xi1 = x1;
                    yi1 = y;
                } else {
                    xi2 = x1;
                    yi2 = y;
                }
                n++;
            }
            return n;
        }
        
        // Oblique segment

        // y = p x + q
        double p = (y2-y1)/(x2-x1);
        double q = y1 - p * x1;
        double denom = a + c*p;
        if (denom == 0) {
            return 0;
        }
        double x;
        int n = 0;
        x = (-c*q-d-w)/denom;
        if ((x1 < x) && (x < x2)) {
            xi1 = x;
            yi1 = p*x + q;
            n = 1;
        }
        x += 2*w/denom;
        if ((x1 < x) && (x < x2)) {
            double y = p*x + q;
            if (n == 0) {
                xi1 = x;
                yi1 = y;
            } else {
                xi2 = x;
                yi2 = y;
            }
            n++;
        }
        return n;
    }

    // ELLIPSE

    // Vertical segment

    if (x1 == x2) {
        // A y^2 + 2 B y + C = 0
        double A = b; // > 0
        double B = c*x1 + e;
        double C = a*x1*x1 + 2*d*x1 + f;
        double discr = B*B - A*C;
        if (discr <= 0) return 0;
        double sqrt_discr = sqrt(discr);
        if (y2 < y1) {
            double tmp = y1; y1 = y2; y2 = tmp;
        }
        int n = 0;
        double y = (-B - sqrt_discr)/A;
        if ((y1 < y) && (y < y2)) {
            xi1 = x1;
            yi1 = y;
            n = 1;
        }
        y = (-B + sqrt_discr)/A;
        if ((y1 < y) && (y < y2)) {
            if (n == 0) {
                xi1 = x1;
                yi1 = y;
            } else {
                xi2 = x1;
                yi2 = y;
            }
            n++;
        }
        return n;
    }
    
    // Oblique segment

    // y = p x + q
    double p = (y2-y1)/(x2-x1);
    double q = y1 - p * x1;
    // A x^2 + 2 B x + C = 0
    // A = (sqrt(b) p + c/sqrt(b))^2 + (ab-c^2)/b > 0
    double A = a + b*p*p + 2*c*p;
    double B = b*p*q + c*q + d + e*p;
    double C = b*q*q + 2*e*q + f;
    double discr = B*B - A*C;
    if (discr <= 0) return 0;
    double sqrt_discr = sqrt(discr);
    int n = 0;
    double x = (-B - sqrt_discr)/A;
    if ((x1 < x) && (x < x2)) {
        xi1 = x;
        yi1 = p*x + q;
        n = 1;
    }
    x = (-B + sqrt_discr)/A;
    if ((x1 < x) && (x < x2)) {
        double y = p*x + q;
        if (n == 0) {
            xi1 = x;
            yi1 = y;
        } else {
            xi2 = x;
            yi2 = y;
        }
        n++;
    }
    return n;
}

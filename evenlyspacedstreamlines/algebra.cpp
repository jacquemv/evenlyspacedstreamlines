#include <math.h>
#include <stdio.h>

//-----------------------------------------------------------------------------
inline void vcopy(double* y, double* x)
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
}

//-----------------------------------------------------------------------------
inline void vswap(double* y, double* x)
{
    double tmp;
    #define SWAP(a, b) tmp = a; a = b; b = tmp;
    SWAP(x[0], y[0]); SWAP(x[1], y[1]); SWAP(x[2], y[2]); 
    #undef SWAP
}

//-----------------------------------------------------------------------------
inline void vcross(double* y, double* r, double* s)
{
    y[0] = r[1]*s[2]-r[2]*s[1];
    y[1] = r[2]*s[0]-r[0]*s[2];
    y[2] = r[0]*s[1]-r[1]*s[0];
}

//-----------------------------------------------------------------------------
inline double vdot(double* x, double* y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//-----------------------------------------------------------------------------
inline double vnorm2(double* x)
{
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

//-----------------------------------------------------------------------------
inline double vdist(double* x, double* y)
{
    return sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) 
              + (x[2]-y[2])*(x[2]-y[2]));
}

//-----------------------------------------------------------------------------
inline double vdist2(double* x, double* y)
{
    return ((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) 
              + (x[2]-y[2])*(x[2]-y[2]));
}

//-----------------------------------------------------------------------------
inline void vscale(double* x, double s)
{
    x[0] *= s;
    x[1] *= s;
    x[2] *= s;
}

//-----------------------------------------------------------------------------
inline void vadd(double* x, double s, double* y)
{
    x[0] += s*y[0];
    x[1] += s*y[1];
    x[2] += s*y[2];
}

//-----------------------------------------------------------------------------
inline void vcopyadd(double* z, double* x, double s, double* y)
{
    z[0] = x[0] + s*y[0];
    z[1] = x[1] + s*y[1];
    z[2] = x[2] + s*y[2];
}

//-----------------------------------------------------------------------------
inline void vlincomb(double* z, double a, double* x, double b, double *y)
{
    z[0] = a*x[0] + b*y[0];
    z[1] = a*x[1] + b*y[1];
    z[2] = a*x[2] + b*y[2];
}

//-----------------------------------------------------------------------------
inline void vdiff(double* z, double* x, double* y)
{
    z[0] = y[0] - x[0];
    z[1] = y[1] - x[1];
    z[2] = y[2] - x[2];
}

//-----------------------------------------------------------------------------
inline void vnormalize(double* x)
{
    vscale(x, 1./sqrt(vnorm2(x)));
}

//-----------------------------------------------------------------------------
inline void vrmcomp(double* x, double* n)
{
    vadd(x, -vdot(x, n)/vnorm2(n), n);
}

//-----------------------------------------------------------------------------
void vprint(double *x)
{
    printf("%.4f %.4f %.4f\n", x[0], x[1], x[2]);
}

//-----------------------------------------------------------------------------
void solve2x2(double a11, double a12, double a21, double a22, 
              double b1, double b2, double& x1, double& x2)
{
    double delta = a11*a22-a12*a21;
    if (delta == 0) return;
    x1 = (b1*a22 - b2*a12)/delta;
    x2 = (b2*a11 - b1*a21)/delta;
}
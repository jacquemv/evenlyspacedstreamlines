#ifndef ALGEBRA_H_
#define ALGEBRA_H_

// basic vector operations
inline void vcopy(double* y, double* x);
inline void vswap(double* y, double* x);
inline void vcross(double* y, double* r, double* s);
inline double vdot(double* x, double* y);
inline double vnorm2(double* x);
inline double vdist(double* x, double* y);
inline double vdist2(double* x, double* y);
inline void vscale(double* x, double s);
inline void vadd(double* x, double s, double* y);
inline void vcopyadd(double* z, double* x, double s, double* y);
inline void vlincomb(double* z, double a, double* x, double b, double *y);
inline void vdiff(double* z, double* x, double* y);
inline void vnormalize(double* x);
inline void vrmcomp(double* x, double* n);
void solve2x2(double a11, double a12, double a21, double a22, 
              double b1, double b2, double& x1, double& x2);
void vprint(double *x);

#endif
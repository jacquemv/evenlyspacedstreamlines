#include "algebra.h"

//-----------------------------------------------------------------------------
double segment_distance2(double *p0, double *p1, double *q0, double *q1)
// distance between the segments P0-P1 and Q0-Q1
// Robust Computation of Distance Between Line Segments, David Eberly, 2020
{
	double a, b, c, d, e, f, det, bte, ctd, ate, btd, s, t;

	// a = (P1 - P0) . (P1 - P0)
	a = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]); 
	// b = (P1 - P0) . (Q1 - Q0)
	b = (p1[0]-p0[0])*(q1[0]-q0[0]) + (p1[1]-p0[1])*(q1[1]-q0[1]) + (p1[2]-p0[2])*(q1[2]-q0[2]); 
	// c = (Q1 - Q0) . (Q1 - Q0)
	c = (q1[0]-q0[0])*(q1[0]-q0[0]) + (q1[1]-q0[1])*(q1[1]-q0[1]) + (q1[2]-q0[2])*(q1[2]-q0[2]); 
	// d = (P1 - P0) . (P0 - Q0)
	d = (p1[0]-p0[0])*(p0[0]-q0[0]) + (p1[1]-p0[1])*(p0[1]-q0[1]) + (p1[2]-p0[2])*(p0[2]-q0[2]); 
	// e = (Q1 - Q0) . (P0 - Q0)
	e = (q1[0]-q0[0])*(p0[0]-q0[0]) + (q1[1]-q0[1])*(p0[1]-q0[1]) + (q1[2]-q0[2])*(p0[2]-q0[2]); 
	// f = (P0 - Q0) . (P0 - Q0)
	f = (p0[0]-q0[0])*(p0[0]-q0[0]) + (p0[1]-q0[1])*(p0[1]-q0[1]) + (p0[2]-q0[2])*(p0[2]-q0[2]); 

	det=a*c-b*b;
	if (det>0) { //non parallel segments
		bte=b*e,ctd=c*d;
		if (bte<=ctd) { //s<=0
			if(e<=0) { //t<=0(region6)
				s=(-d>=a?1:(-d>0?-d/a:0));t=0;
			}
			else if (e<c) { //0<t<1(region5)
				s=0;t=e/c;
			} else { //t>=1(region4)
				s=(b-d>=a?1:(b-d>0?(b-d)/a:0));t=1;
			}
		} else { //s>0
			s=bte-ctd;
			if (s>=det) { //s>=1
				if (b+e<=0) { //t<=0(region8)
					s=(-d<=0?0:(-d<a?-d/a:1));t=0;
				} else if (b+e<c) {//0<t<1(region1)
					s=1;t=(b+e)/c;
				} else { //t>=1(region2)
					s=(b-d<=0?0:(b-d<a?(b-d)/a:1));t=1;
				}
			} else { //0<s<1
				ate=a*e,btd=b*d;
				if (ate<=btd) { //t<=0(region7)
					s=(-d<=0?0:(-d>=a?1:-d/a));t=0;
				} else {//t>0
					t=ate-btd;
					if (t>=det) { //t>=1(region3)
						s=(b-d<=0?0:(b-d>=a?1:(b-d)/a));t=1;
					} else { //0<t<1(region0)
						s/=det;
						t/=det;
					}
				}
			}
		}
	} else { //parallel segments
		if (e<=0) {
			s=(-d<=0?0:(-d>=a?1:-d/a));t=0;
		} else if (e>=c) {
			s=(b-d<=0?0:(b-d>=a?1:(b-d)/a));t=1;
		} else {
			s=0;t=e/c;
		}
	}

	double dist2 = a*s*s - 2*b*s*t + c*t*t + 2*d*s - 2*e*t + f;
	if (dist2 < 0) dist2 = 0;
	return dist2;
}

//-----------------------------------------------------------------------------
bool segments_close_enough(double *p0, double *p1, double *q0, double *q1, 
						   double radius2)
{
	double a, b, c, d, e, f, det, bte, ctd, ate, btd, s, t;

	// f = (P0 - Q0) . (P0 - Q0)
	f = (p0[0]-q0[0])*(p0[0]-q0[0]) + (p0[1]-q0[1])*(p0[1]-q0[1]) + (p0[2]-q0[2])*(p0[2]-q0[2]);
	if (f < radius2) return 1;

	// c = (Q1 - Q0) . (Q1 - Q0)
	c = (q1[0]-q0[0])*(q1[0]-q0[0]) + (q1[1]-q0[1])*(q1[1]-q0[1]) + (q1[2]-q0[2])*(q1[2]-q0[2]); 
	// e = (Q1 - Q0) . (P0 - Q0)
	e = (q1[0]-q0[0])*(p0[0]-q0[0]) + (q1[1]-q0[1])*(p0[1]-q0[1]) + (q1[2]-q0[2])*(p0[2]-q0[2]); 
	if (c-2*e+f < radius2) return 1;

	// a = (P1 - P0) . (P1 - P0)
	a = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
	// d = (P1 - P0) . (P0 - Q0)
	d = (p1[0]-p0[0])*(p0[0]-q0[0]) + (p1[1]-p0[1])*(p0[1]-q0[1]) + (p1[2]-p0[2])*(p0[2]-q0[2]);
	if (a+2*d+f < radius2) return 1;

	// b = (P1 - P0) . (Q1 - Q0)
	b = (p1[0]-p0[0])*(q1[0]-q0[0]) + (p1[1]-p0[1])*(q1[1]-q0[1]) + (p1[2]-p0[2])*(q1[2]-q0[2]);
	if (c - 2*e + a + 2*d - 2*b + f < radius2) return 1;
	
	det=a*c-b*b;
	if (det>0) { //non parallel segments
		bte=b*e,ctd=c*d;
		if (bte<=ctd) { //s<=0
			if(e<=0) { //t<=0(region6)
				s=(-d>=a?1:(-d>0?-d/a:0));t=0;
			}
			else if (e<c) { //0<t<1(region5)
				s=0;t=e/c;
			} else { //t>=1(region4)
				s=(b-d>=a?1:(b-d>0?(b-d)/a:0));t=1;
			}
		} else { //s>0
			s=bte-ctd;
			if (s>=det) { //s>=1
				if (b+e<=0) { //t<=0(region8)
					s=(-d<=0?0:(-d<a?-d/a:1));t=0;
				} else if (b+e<c) {//0<t<1(region1)
					s=1;t=(b+e)/c;
				} else { //t>=1(region2)
					s=(b-d<=0?0:(b-d<a?(b-d)/a:1));t=1;
				}
			} else { //0<s<1
				ate=a*e,btd=b*d;
				if (ate<=btd) { //t<=0(region7)
					s=(-d<=0?0:(-d>=a?1:-d/a));t=0;
				} else {//t>0
					t=ate-btd;
					if (t>=det) { //t>=1(region3)
						s=(b-d<=0?0:(b-d>=a?1:(b-d)/a));t=1;
					} else { //0<t<1(region0)
						s/=det;
						t/=det;
					}
				}
			}
		}
	} else { //parallel segments
		if (e<=0) {
			s=(-d<=0?0:(-d>=a?1:-d/a));t=0;
		} else if (e>=c) {
			s=(b-d<=0?0:(b-d>=a?1:(b-d)/a));t=1;
		} else {
			s=0;t=e/c;
		}
	}

	double dist2 = a*s*s - 2*b*s*t + c*t*t + 2*d*s - 2*e*t + f;
	return dist2 < radius2;
}

//-----------------------------------------------------------------------------
bool intersect_sphere(double* p, double* q, double* a, double r, 
					  double& tmin, double& tmax)
// intersection between segment PQ and the sphere (a, r)
// the result is the interval [tmin,tmax] inter [0,1]
// returns 0 if there is no intersection
{
	double u[3] = {q[0]-p[0], q[1]-p[1], q[2]-p[2]};
	double A = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
	double B = u[0]*(p[0]-a[0]) + u[1]*(p[1]-a[1]) + u[2]*(p[2]-a[2]);
	double C = (p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]) 
			+ (p[2]-a[2])*(p[2]-a[2]) - r*r;
	double D2 = B*B - A*C;
	if (D2 <= 0 || A == 0) return 0;
	double D = sqrt(D2);
	tmin = (-B-D)/A;
	tmax = (-B+D)/A;
	if ((tmin >= tmax) || (tmin >= 1) || (tmax <= 0)) return 0;
	return 1;
}

//-----------------------------------------------------------------------------
bool intersect_cylinder(double* p, double* q, double* a, double* b, double r, 
						double& tmin, double& tmax)
// intersection between segment PQ and the cylinder of axis (a, b) and radius r
// the result is the interval [tmin,tmax] inter [0,1]
// returns 0 if there is no intersection
{
	double u[3] = {q[0]-p[0], q[1]-p[1], q[2]-p[2]};
	double v[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	if (v2 == 0) return 0;
	double uv_v2 = (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])/v2;
	double w[3] = {uv_v2*v[0]-u[0], uv_v2*v[1]-u[1], uv_v2*v[2]-u[2]};
	double m = ((p[0]-a[0])*v[0] + (p[1]-a[1])*v[1] + (p[2]-a[2])*v[2])/v2;
	double proj[3] = {a[0]-p[0]+m*v[0], a[1]-p[1]+m*v[1], a[2]-p[2]+m*v[2]};
	double C = proj[0]*proj[0] + proj[1]*proj[1] + proj[2]*proj[2] - r*r;
	double B = w[0]*proj[0] + w[1]*proj[1] + w[2]*proj[2];
	double A = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
	double D2 = B*B - A*C;
	if (D2 < 0) return 0;
	double t1, t2, t3, t4;
	if (A < 1e-12) { // so B = D = 0
		if (C >= 0) return 0;
		t1 = 0;
		t2 = 1;
	} else {
		double D = sqrt(D2);
		t1 = (-B-D)/A;
		t2 = (-B+D)/A;
	}
	if (uv_v2 == 0) {
		if ((m <= 0) || (m >= 1)) return 0;
		t3 = 0;
		t4 = 1;
	} else {
		t3 = -m/uv_v2;
		t4 = (1-m)/uv_v2;
		if (t4 < t3) {
			double tmp = t3;
			t3 = t4;
			t4 = tmp;
		}
	}
	tmin = t1 > t3 ? t1 : t3;
	tmax = t2 > t4 ? t4 : t2;
	if ((tmin >= tmax) || (tmin >= 1) || (tmax <= 0)) return 0;
	return 1;
}

//-----------------------------------------------------------------------------
int connect_segments(double* A, double* B, double* C, double* D)
// AC closest -> 0x00, AD closest -> 0x01, 
// BC closest -> 0x10, BD closest -> 0x11
{
    double d2[4] = {vdist2(A, C), vdist2(A, D), vdist2(B, C), vdist2(B, D)};
    double d2min = d2[0];
    int code = 0x00;
    if (d2[1] < d2min) { code = 0x01; d2min = d2[1]; }
    if (d2[2] < d2min) { code = 0x10; d2min = d2[2]; }
    if (d2[3] < d2min) { code = 0x11; }
    return code;
}
#ifndef GEOMETRY_H_
#define GEOMETRY_H_

// basic geometry routines
double segment_distance2(double *p0, double *p1, double *q0, double *q1);
bool segments_close_enough(double *p0, double *p1, double *q0, double *q1, 
						   double radius2);
bool intersect_sphere(double* p, double* q, double* a, double r, 
					  double& tmin, double& tmax);
bool intersect_cylinder(double* p, double* q, double* a, double* b, double r, 
						double& tmin, double& tmax);
int connect_segments(double* A, double* B, double* C, double* D);

#endif
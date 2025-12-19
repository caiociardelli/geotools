/*
 GeoTools

 Author: Caio Ciardelli, Northwestern University, October 2025

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

----------------------------------------------------------------------------------------------- */

#include <math.h>
#include "constants.h"
#include "exmath.h"
#include "coordinates.h"

double angle (struct Triplet *p, struct Triplet *q)
{
  /* Computes the angle between two 3D vectors in radians. */
  normalize (p);
  normalize (q);
  
  return acos (dot (p, q));
}

void rThetaPhi2XYZ (double r, double theta, double phi,
                    double *x, double *y, double *z)
{
  /* Converts spherical coordinates to Cartesian */
  *x = r * sin (theta) * cos (phi);
  *y = r * sin (theta) * sin (phi);
  *z = r * cos (theta);
}

void xYZ2RThetaPhi (double x, double y, double z,
                    double *r, double *theta, double *phi)
{
  /* Converts Cartesian coordinates to spherical */
  *r = sqrt (square (x) + square (y) + square (z));

  if (z > -ANGLE_TOLERANCE && z <= 0.0) z = -ANGLE_TOLERANCE;
  if (z <  ANGLE_TOLERANCE && z >= 0.0) z =  ANGLE_TOLERANCE;

  *theta = atan2 (sqrt (square (x) + square (y)), z);

  if (x > -ANGLE_TOLERANCE && x <= 0.0) x = -ANGLE_TOLERANCE;
  if (x <  ANGLE_TOLERANCE && x >= 0.0) x =  ANGLE_TOLERANCE;

  *phi = atan2 (y, x);
}

void capCoordinates (double theta, double phi,
                     double *lat, double *lon)
{
  /* Caps and converts spherical coordinates (theta, phi)
     to latitude and longitude within bounds. */
  *lat = 90.0 - rad2Degree (theta);
  *lon = rad2Degree (phi);

  if      (*lat > LAT_MAX) *lat = LAT_MAX;
  else if (*lat < LAT_MIN) *lat = LAT_MIN;
  if      (*lon > LON_MAX) *lon = LON_MAX;
  else if (*lon < LON_MIN) *lon = LON_MIN;
}

struct Triplet midPoint (struct Triplet *p1, struct Triplet *p2)
{
  /* Computes the midpoint of two 3D points represented by Triplet structures. */
  struct Triplet pm;

  pm.x = 0.5 * (p1->x + p2->x);
  pm.y = 0.5 * (p1->y + p2->y);
  pm.z = 0.5 * (p1->z + p2->z);

  return pm;
}

struct Triplet normalVector (struct Triplet *p1,
                             struct Triplet *p2,
                             struct Triplet *p3)
{
  /* Computes the normalized normal vector of a triangular facet. */
  struct Triplet u, v;

  u.x = p3->x - p1->x; u.y = p3->y - p1->y; u.z = p3->z - p1->z;
  v.x = p2->x - p1->x; v.y = p2->y - p1->y; v.z = p2->z - p1->z;

  struct Triplet n = cross (&u, &v);

  normalize (&n);

  return n;
}

struct Triplet project (struct Triplet *p,
                        struct Triplet *p1,
                        struct Triplet *p2,
                        struct Triplet *p3)
{
  /* Projects a point onto the plane defined by three points. */
  struct Triplet n = normalVector (p1, p2, p3);

  double a = n.x;
  double b = n.y;
  double c = n.z;
  double d = -(a * p1->x + b * p1->y + c * p1->z);

  double k = fabs (-d / (a * p->x + b * p->y + c * p->z + EPSILON));

  struct Triplet q;

  q.x = k * p->x;
  q.y = k * p->y;
  q.z = k * p->z;

  return q;
}

void baricentric (struct Triplet *p,
                  struct Triplet *p1,
                  struct Triplet *p2,
                  struct Triplet *p3,
                  double *u, double *v, double *w)
{
  /* Computes barycentric coordinates of a point relative to a triangle. */
  double x1 = p1->x, y1 = p1->y, z1 = p1->z;
  double x2 = p2->x, y2 = p2->y, z2 = p2->z;
  double x3 = p3->x, y3 = p3->y, z3 = p3->z;

  /* Compute normal vector of the triangle */
  struct Triplet n = normalVector (p1, p2, p3);

  /* Vectors from triangle vertices to each other and to the point */
  struct Triplet v21, v31, v32;

  v21.x = x2 - x1; v21.y = y2 - y1; v21.z = z2 - z1;
  v31.x = x3 - x1; v31.y = y3 - y1; v31.z = z3 - z1;
  v32.x = x3 - x2; v32.y = y3 - y2; v32.z = z3 - z2;

  /* Vectors from triangle vertices to the point */
  struct Triplet vp1, vp2;

  vp1.x = p->x - x1; vp1.y = p->y - y1; vp1.z = p->z - z1;
  vp2.x = p->x - x2; vp2.y = p->y - y2; vp2.z = p->z - z2;

  /* Cross products for area calculations */
  struct Triplet w1 = cross (&v21, &v31);
  struct Triplet w2 = cross (&v21, &vp1);
  struct Triplet w3 = cross (&v31, &vp1);
  struct Triplet w4 = cross (&v32, &vp2);

  /* Compute areas for barycentric coordinates */
  double ABC = fabs (dot (&n, &w1));
  double ABP = fabs (dot (&n, &w2));
  double ACP = fabs (dot (&n, &w3));
  double BCP = fabs (dot (&n, &w4));

  /* Inverse of total triangle area for normalization */
  double ABCi = 1.0 / (ABC + EPSILON);

  /* Calculate barycentric coordinates */
  *u = BCP * ABCi; *v = ACP * ABCi; *w = ABP * ABCi;
}

void weights (struct Triplet *p,
              double r1, double r2, double r3,
              double r4, double r5, double r6,
              struct Facet *fct,
              int nv, struct Vertex vtx[nv],
              double *w1, double *w2, double *w3,
              double *w4, double *w5, double *w6)
{
  /* Computes interpolation weights for a point within a prismatic element. */
  struct Triplet p1 = vtx[fct->i1].p;
  struct Triplet p2 = vtx[fct->i2].p;
  struct Triplet p3 = vtx[fct->i3].p;
  struct Triplet p4 = vtx[fct->i1].p;
  struct Triplet p5 = vtx[fct->i2].p;
  struct Triplet p6 = vtx[fct->i3].p;

  /* Scale vertices by radial distances for top and bottom facets */
  p1.x *= r1; p1.y *= r1; p1.z *= r1;
  p2.x *= r2; p2.y *= r2; p2.z *= r2;
  p3.x *= r3; p3.y *= r3; p3.z *= r3;
  p4.x *= r4; p4.y *= r4; p4.z *= r4;
  p5.x *= r5; p5.y *= r5; p5.z *= r5;
  p6.x *= r6; p6.y *= r6; p6.z *= r6;

  /* Project point onto top and bottom facets */
  struct Triplet q = project (p, &p1, &p2, &p3);
  struct Triplet r = project (p, &p4, &p5, &p6);

  double uq, vq, wq;

  /* Compute barycentric coordinates for top facet projection */
  baricentric (&q, &p1, &p2, &p3, &uq, &vq, &wq);

  double ur, vr, wr;

  /* Use top facet coordinates if bottom projection is near origin */
  if (norm (&r) < EPSILON)
  {
    ur = uq; vr = vq; wr = wq;
  }

  /* Compute barycentric coordinates for bottom facet projection */
  else baricentric (&r, &p4, &p5, &p6, &ur, &vr, &wr);

  /* Compute norms for interpolation weighting */
  double nmq = norm (&q);
  double nmp = norm (p);
  double nmr = norm (&r);

  /* Calculate radial interpolation factors */
  double npr = nmp - nmr;
  double nqr = nmq - nmr;

  /* Avoid division by zero */
  if (fabs (nqr) < EPSILON) nqr = EPSILON;

  /* Compute weights for linear interpolation between facets */
  double rq = npr / nqr;
  double rr = 1.0 - rq;

  /* Assign final interpolation weights */
  *w1 = rq * uq; *w2 = rq * vq; *w3 = rq * wq;
  *w4 = rr * ur; *w5 = rr * vr; *w6 = rr * wr;
}
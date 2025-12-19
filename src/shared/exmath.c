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
#include "structs.h"
#include "exmath.h"

int findIntervalIndex (int n, double t[n], double s)
{
  /* Finds the index of the interval containing a given value using binary search. */
  int lo = 0, hi = n - 1;

  while (hi - lo > 1)
  {
    int mid = lo + (hi - lo) / 2;

    if (t[mid] <= s) lo = mid; else hi = mid;
  }

  return lo;
}

bool checkRedundancy (struct Triplet *p, struct Triplet *q)
{
  /* Checks if two 3D points are nearly identical within a tolerance. */
  double sum = fabs (p->x - q->x)
             + fabs (p->y - q->y)
             + fabs (p->z - q->z);

  if (sum < TOLERANCE) return true;

  return false;
}

struct Triplet cross (struct Triplet *u, struct Triplet *v)
{
  /* Computes the cross product of two 3D vectors. */
  struct Triplet w;

  w.x = u->y * v->z - u->z * v->y;
  w.y = u->z * v->x - u->x * v->z;
  w.z = u->x * v->y - u->y * v->x;

  return w;
}

void rotate (int nv, struct Triplet Icv[nv],
             double alpha, double beta, double gamma)
{
  /* Rotates a set of 3D points around the x, y, and z axes by given angles. */
  double Rx[N_DIM][N_DIM];
  double Ry[N_DIM][N_DIM];
  double Rz[N_DIM][N_DIM];

  Rx[0][0] =  1.0;         Rx[0][1] =  0.0;         Rx[0][2] =  0.0;
  Rx[1][0] =  0.0;         Rx[1][1] =  cos (alpha); Rx[1][2] = -sin (alpha);
  Rx[2][0] =  0.0;         Rx[2][1] =  sin (alpha); Rx[2][2] =  cos (alpha);

  Ry[0][0] =  cos (beta);  Ry[0][1] =  0.0;         Ry[0][2] =  sin (beta);
  Ry[1][0] =  0.0;         Ry[1][1] =  1.0;         Ry[1][2] =  0.0;
  Ry[2][0] = -sin (beta);  Ry[2][1] =  0.0;         Ry[2][2] =  cos (beta);

  Rz[0][0] =  cos (gamma); Rz[0][1] = -sin (gamma); Rz[0][2] =  0.0;
  Rz[1][0] =  sin (gamma); Rz[1][1] =  cos (gamma); Rz[1][2] =  0.0;
  Rz[2][0] =  0.0;         Rz[2][1] =  0.0;         Rz[2][2] =  1.0;

  /* Apply rotation matrices sequentially to each vertex */
  for (int i = 0; i < nv; i++)
  {
    double x = Icv[i].x, y = Icv[i].y, z = Icv[i].z;

    Icv[i].x = Rx[0][0] * x + Rx[0][1] * y + Rx[0][2] * z;
    Icv[i].y = Rx[1][0] * x + Rx[1][1] * y + Rx[1][2] * z;
    Icv[i].z = Rx[2][0] * x + Rx[2][1] * y + Rx[2][2] * z;

    x = Icv[i].x; y = Icv[i].y; z = Icv[i].z;

    Icv[i].x = Ry[0][0] * x + Ry[0][1] * y + Ry[0][2] * z;
    Icv[i].y = Ry[1][0] * x + Ry[1][1] * y + Ry[1][2] * z;
    Icv[i].z = Ry[2][0] * x + Ry[2][1] * y + Ry[2][2] * z;

    x = Icv[i].x; y = Icv[i].y; z = Icv[i].z;

    Icv[i].x = Rz[0][0] * x + Rz[0][1] * y + Rz[0][2] * z;
    Icv[i].y = Rz[1][0] * x + Rz[1][1] * y + Rz[1][2] * z;
    Icv[i].z = Rz[2][0] * x + Rz[2][1] * y + Rz[2][2] * z;
  }
}

void computeDistanceVector (struct Triplet *p,
                            struct Triplet *q,
                            struct Triplet *r)
{
  /* Computes the distance vector between two points and normalizes it. */
  /* Compute the differences in x, y, and z coordinates */
  r->x = q->x - p->x;
  r->y = q->y - p->y;
  r->z = q->z - p->z;

  /* Calculate the squared magnitude of the vector */
  double r2 = square (r->x) +
              square (r->y) +
              square (r->z);

  /* Compute the normalization factor with a small epsilon for stability */
  double nrm_factor = sqrt (r2) * r2 + EPSILON;

  /* Normalize the x, y, and z components of the vector */
  r->x /= nrm_factor;
  r->y /= nrm_factor;
  r->z /= nrm_factor;
}

void createNormalVectors (int np, int nt,
                          struct Triplet Tn[np][nt])
{
  /* Creates normal vectors for a grid of points. */
  double r = -1;

  double dp = 2 * PI / (np - 1);
  double dt = PI / (nt - 1);

  double phi = -PI;

  for (int m = 0; m < np; m++)
  {
    double theta = 0;

    for (int n = 0; n < nt; n++)
    {
      /* Convert spherical coordinates to Cartesian for normal vector */
      rThetaPhi2XYZ (r, theta, phi,
                     &Tn[m][n].x,
                     &Tn[m][n].y,
                     &Tn[m][n].z);

      theta += dt;
    }

    phi += dp;
  }
}

double tripletDistance (struct Triplet *a, struct Triplet *b)
{
  /* Computes the Euclidean distance between two 3D points. */
  double dx = a->x - b->x;
  double dy = a->y - b->y;
  double dz = a->z - b->z;

  return sqrt (square (dx) + square (dy) + square (dz));
}

double volumeTetrahedron (struct Tetrahedron *th)
{
  /* Computes the volume of a tetrahedron using the scalar triple product. */
  struct Triplet v1, v2, v3;

  /* Calculate vector from p1 to p2 */
  v1.x = th->p2.x - th->p1.x;
  v1.y = th->p2.y - th->p1.y;
  v1.z = th->p2.z - th->p1.z;

  /* Calculate vector from p1 to p3 */
  v2.x = th->p3.x - th->p1.x;
  v2.y = th->p3.y - th->p1.y;
  v2.z = th->p3.z - th->p1.z;

  /* Calculate vector from p1 to p4 */
  v3.x = th->p4.x - th->p1.x;
  v3.y = th->p4.y - th->p1.y;
  v3.z = th->p4.z - th->p1.z;

  /* Compute the cross product of v1 and v2 */
  struct Triplet v12 = cross (&v1, &v2);

  /* Return the absolute volume divided by 6 */
  return fabs (dot (&v12, &v3)) / 6;
}
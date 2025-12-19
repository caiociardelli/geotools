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

#ifndef EXMATH_H
#define EXMATH_H

#include "constants.h"
#include "structs.h"
#include "coordinates.h"

static inline double square (double v)
{
  /* Computes the square of a scalar value. */
  return v * v;
}

static inline double norm (struct Triplet *p)
{
  /* Computes the Euclidean norm of a 3D vector. */
  return sqrt (p->x * p->x + p->y * p->y + p->z * p->z);
}

static inline double goldenRatio (void)
{
  /* Returns the golden ratio, (1 + sqrt(5)) / 2. */
  return (1.0 + sqrt (5)) / 2;
}

static inline void normalize (struct Triplet *p)
{
  /* Normalizes a 3D vector to unit length. */
  double inorm = 1.0 / (norm (p) + EPSILON);

  p->x *= inorm;
  p->y *= inorm;
  p->z *= inorm;
}

static inline double dot (struct Triplet *u, struct Triplet *v)
{
  /* Computes the dot product of two 3D vectors. */
  return u->x * v->x + u->y * v->y + u->z * v->z;
}

int findIntervalIndex (int n, double t[n], double s);

bool checkRedundancy (struct Triplet *p, struct Triplet *q);

struct Triplet cross (struct Triplet *u, struct Triplet *v);

void rotate (int nv, struct Triplet Icv[nv],
             double alpha, double beta, double gamma);
void computeDistanceVector (struct Triplet *p,
                            struct Triplet *q,
                            struct Triplet *r);
void createNormalVectors (int np, int nt,
                          struct Triplet Tn[np][nt]);

double tripletDistance (struct Triplet *a, struct Triplet *b);
double volumeTetrahedron (struct Tetrahedron *th);
#endif
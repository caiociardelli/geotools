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

#ifndef COORDINATES_H
#define COORDINATES_H

#include "constants.h"
#include "structs.h"

static inline double rad2Degree (double v)
{
  /* Converts radians to degrees. */
  return v * TO_DEGREE;
}

static inline double degree2Rad (double v)
{
  /* Converts degrees to radians. */
  return v * TO_RADIAN;
}

double angle (struct Triplet *p, struct Triplet *q);

void rThetaPhi2XYZ (double r, double theta, double phi,
                    double *x, double *y, double *z);
void xYZ2RThetaPhi (double x, double y, double z,
                    double *r, double *theta, double *phi);
void capCoordinates (double theta, double phi,
                     double *lat, double *lon);

struct Triplet midPoint (struct Triplet *p1, struct Triplet *p2);
struct Triplet normalVector (struct Triplet *p1,
                             struct Triplet *p2,
                             struct Triplet *p3);
struct Triplet project (struct Triplet *p,
                        struct Triplet *p1,
                        struct Triplet *p2,
                        struct Triplet *p3);
void baricentric (struct Triplet *p,
                  struct Triplet *p1,
                  struct Triplet *p2,
                  struct Triplet *p3,
                  double *u, double *v, double *w);
void weights (struct Triplet *p,
              double r1, double r2, double r3,
              double r4, double r5, double r6,
              struct Facet *fct,
              int nv, struct Vertex vtx[nv],
              double *w1, double *w2, double *w3,
              double *w4, double *w5, double *w6);
#endif
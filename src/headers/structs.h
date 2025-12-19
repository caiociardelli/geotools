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

#ifndef STRUCTS_H
#define STRUCTS_H

#include "constants.h"

struct Indices
{
  /* Stores indices for a triangular facet. */
  short i1;
  short i2;
  short i3;
};

struct Triplet
{
  /* Represents a 3D point with x, y, z coordinates. */
  double x;
  double y;
  double z;
};

struct FacetPoints
{
  /* Defines a facet with points, neighbors, and spherical normal. */
  int nnf;                         /* Number of neighboring facets */
  int nbf[N_MAX_NGB_FACETS_MACRO]; /* Array of neighboring facet indices */

  struct Triplet p1;               /* First vertex of the facet */
  struct Triplet p2;               /* Second vertex of the facet */
  struct Triplet p3;               /* Third vertex of the facet */

  struct Triplet n_sph;            /* Spherical normal vector of the facet */

  double min_dot;                  /* Minimum dot product for facet projection */
};

struct IndexGroup
{
  /* Holds indices and IDs for neighboring triangles. */
  int id[N_MAX_NGB_TRIAG_MACRO]; /* Array of neighbor triangle IDs */
  int ii[N_MAX_NGB_TRIAG_MACRO]; /* Array of neighbor triangle indices */
};

struct Vertex
{
  /* Represents a vertex with position and neighbor information. */
  int nnv;                           /* Number of neighboring vertices */
  int nbv[N_MAX_NGB_VERTICES_MACRO]; /* Array of neighboring vertex indices */

  struct Triplet p;                  /* 3D coordinates of the vertex */
};

struct Facet
{
  /* Defines a mesh facet with connectivity and spherical normal. */
  int nnf;                         /* Number of neighboring facets */
  int nbf[N_MAX_NGB_FACETS_MACRO]; /* Array of neighboring facet indices */

  int nbf_12;                      /* Neighbor facet sharing vertices i1 and i2 */
  int nbf_13;                      /* Neighbor facet sharing vertices i1 and i3 */
  int nbf_23;                      /* Neighbor facet sharing vertices i2 and i3 */

  int i1;                          /* Index of first vertex */
  int i2;                          /* Index of second vertex */
  int i3;                          /* Index of third vertex */

  struct Triplet n_sph;            /* Spherical normal vector of the facet */

  double min_dot;                  /* Minimum dot product for facet projection */
};

struct Element
{
  /* Represents a prismatic element with densities and coordinates at six vertices. */
  double rho1; /* Density at vertex 1 */
  double rho2; /* Density at vertex 2 */
  double rho3; /* Density at vertex 3 */
  double rho4; /* Density at vertex 4 */
  double rho5; /* Density at vertex 5 */
  double rho6; /* Density at vertex 6 */

  struct Triplet p1; /* Coordinate of vertex 1 */
  struct Triplet p2; /* Coordinate of vertex 2 */
  struct Triplet p3; /* Coordinate of vertex 3 */
  struct Triplet p4; /* Coordinate of vertex 4 */
  struct Triplet p5; /* Coordinate of vertex 5 */
  struct Triplet p6; /* Coordinate of vertex 6 */
};

struct Tetrahedron
{
  /* Represents a tetrahedron with coordinates of four vertices. */
  struct Triplet p1; /* Coordinate of vertex 1 */
  struct Triplet p2; /* Coordinate of vertex 2 */
  struct Triplet p3; /* Coordinate of vertex 3 */
  struct Triplet p4; /* Coordinate of vertex 4 */
};

struct Prism
{
  /* Represents a prism with a skip flag, mass, and center of mass. */
  bool skip; /* Flag indicating whether to skip this prism in calculations */

  double dm; /* Mass of the prism */

  struct Triplet cm; /* Center of mass of the prism */
};
#endif
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

/* Do not change anything in this file unless you know really well what you are doing! */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <stdbool.h>

#define N_MAX_NGB_FACETS_MACRO   13                             /* Maximum number of neighboring facets */
#define N_MAX_NGB_TRIAG_MACRO    6                              /* Maximum number of neighboring triangles */
#define N_MAX_NGB_VERTICES_MACRO 6                              /* Maximum number of neighboring vertices */

static const int MAX_STRING_LEN = 200;                          /* Buffer size for short strings */
static const int MAX_PATH_LEN   = 300;                          /* Buffer size for file paths */

static const int N_CORNERS          = 3;                        /* Number of corners in a triangle */
static const int N_DIM              = 3;                        /* Number of dimensions */
static const int N_MAX_NGB_FACETS   = N_MAX_NGB_FACETS_MACRO;   /* Maximum number of neighboring facets */
static const int N_MAX_NGB_TRIAG    = N_MAX_NGB_TRIAG_MACRO;    /* Maximum number of neighboring triangles */
static const int N_MAX_NGB_VERTICES = N_MAX_NGB_VERTICES_MACRO; /* Maximum number of neighboring vertices */
static const int N_LAT_ETOPO        = 721;                      /* Number of latitude points in ETOPO data */
static const int N_LON_ETOPO        = 1441;                     /* Number of longitude points in ETOPO data */
static const int N_TETRAHEDRA       = 32;                       /* Number of tetrahedra in the mesh */

static const bool INCORPORATE_SURFACE_TOPOGRAPHY = true;        /* Flag to incorporate surface topography */
static const bool INCORPORATE_3D_MODEL           = false;       /* Flag to incorporate 3D velocity and density model */
static const bool DEBUG_MESH                     = true;        /* Flag to print edges and facets for plotting */

static const double EARTH_RADIUS    = 6371.0;                   /* Earth's radius in kilometers */
static const double PI              = 3.14159265358979323846;   /* Pi approximation */
static const double TO_RADIAN       = PI / 180.0;               /* Constant to convert from degrees to radians */
static const double TO_DEGREE       = 180.0 / PI;               /* Constant to convert from radians to degrees */
static const double EPSILON         = 1E-15;                    /* Small value for numerical stability */
static const double THRESHOLD       = 1E-9;                     /* Threshold for radius comparison */
static const double TOLERANCE       = 1E-3;                     /* Tolerance for point redundancy check */
static const double LEEWAY          = 1.00001;                  /* Leeway for barycentric coordinate checks */
static const double DELTA           = 0.1;                      /* Spacing for ray path resampling */
static const double LAT_MAX         = 89.999999;                /* Maximum latitude bound */
static const double LAT_MIN         = -89.999999;               /* Minimum latitude bound */
static const double LON_MAX         =  179.999999;              /* Maximum longitude bound */
static const double LON_MIN         = -179.999999;              /* Minimum longitude bound */
static const double ANGLE_TOLERANCE = 1E-10;                    /* Tolerance for angular comparisons */
static const double MIN_DOT_NODES   = 9.9999999999999E-1;       /* Minimum dot product for node matching */
static const double SKIP_VOLUME     = 1E-7;                     /* Threshold for skipping small volumes */
static const double G               = 6.67384E-17;              /* Gravitational constant in m^3 kg^-1 s^-2 */
#endif
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

-----------------------------------------------------------------------------------------------

 GRAVITY

 USAGE
   mpiexec -n N ./gravity DH HEIGHT

 EXAMPLE
   mpiexec -n 4 ./gravity 1.0 255

 COMMAND-LINE ARGUMENTS
   DH                    - grid spacing in degrees
   HEIGHT                - observation height in kilometers

 DESCRIPTION
   Computes the Earth's gravitational field using a spherical mesh and velocity model from
   'input.model'. Accounts for surface topography and observation height. Reads mesh and model
   data from binary files (facets.bin, vertices.bin, r.bin, rho.bin, vp.bin, vs.bin) and
   outputs gravity components (x, y, z) and normal component to text files (G_x.dat,
   G_y.dat, G_z.dat, G_normal.dat). Uses MPI for parallel computation.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "constants.h"
#include "structs.h"
#include "io.h"
#include "exmath.h"
#include "coordinates.h"
#include "tesselation.h"
#include "progress.h"

static void partitionFacets (int ic, int nc,
                             int nf, int *nf_l,
                             int *shift)
{
  /* Partitions facets among processes for parallel computation. */
  int nf_lr = nf % nc, nf_lv[nc];

  *nf_l = nf / nc; *shift = 0;

  for (int i = 0; i < nc; i++)
  {
    /* Assign number of facets per process with load balancing */
    nf_lv[i] = (i < nf_lr) ? *nf_l + 1 : *nf_l;

    if (i < ic) *shift += nf_lv[i];
  }

  *nf_l = nf_lv[ic];
}

static void initiallizeArray (int np, int nt,
                              struct Triplet Gl[np][nt])
{
  /* Initializes a 2D array of Triplet structures to zero. */
  for (int m = 0; m < np; m++)

    for (int n = 0; n < nt; n++)
    {
      /* Set x-component to zero */
      Gl[m][n].x = 0;
      /* Set y-component to zero */
      Gl[m][n].y = 0;
      /* Set z-component to zero */
      Gl[m][n].z = 0;
    }
}

static void createOutputGrid (double height,
                              int np, int nt,
                              struct Triplet To[np][nt])
{
  /* Creates a spherical grid of observation points at a specified height. */
  double r = EARTH_RADIUS + height;

  double dp = 2 * PI / (np - 1);
  double dt = PI / (nt - 1);

  double phi = -PI;

  for (int m = 0; m < np; m++)
  {
    double theta = 0;

    for (int n = 0; n < nt; n++)
    {
      /* Convert spherical coordinates to Cartesian */
      rThetaPhi2XYZ (r, theta, phi,
                     &To[m][n].x,
                     &To[m][n].y,
                     &To[m][n].z);

      theta += dt;
    }

    phi += dp;
  }
}

static void createMesh (int nv, int ns,
                        double r[nv][ns],
                        double rho[nv][ns],
                        int nf, int nf_l,
                        int shift,
                        struct Facet fct[nf],
                        struct Vertex vtx[nv],
                        int nl,
                        struct Element elm[nl][nf_l])
{
  /* Creates a mesh of prismatic elements from facet and vertex data with scaled
     radii and densities. */
  for (int j = 0; j < nf_l; j++)
  {
    /* Extract vertex indices for the current facet with offset */
    int i1 = fct[j + shift].i1;
    int i2 = fct[j + shift].i2;
    int i3 = fct[j + shift].i3;

    for (int i = 1; i < ns; i++)
    {
      /* Retrieve radius values for the previous and current shell at each vertex */
      double r1 = r[i1][i - 1];
      double r2 = r[i2][i - 1];
      double r3 = r[i3][i - 1];
      double r4 = r[i1][i];
      double r5 = r[i2][i];
      double r6 = r[i3][i];

      /* Set coordinates for vertex 1 of the element using previous shell radius */
      elm[i][j].p1.x = r1 * vtx[i1].p.x;
      elm[i][j].p1.y = r1 * vtx[i1].p.y;
      elm[i][j].p1.z = r1 * vtx[i1].p.z;

      /* Set coordinates for vertex 2 of the element using previous shell radius */
      elm[i][j].p2.x = r2 * vtx[i2].p.x;
      elm[i][j].p2.y = r2 * vtx[i2].p.y;
      elm[i][j].p2.z = r2 * vtx[i2].p.z;

      /* Set coordinates for vertex 3 of the element using previous shell radius */
      elm[i][j].p3.x = r3 * vtx[i3].p.x;
      elm[i][j].p3.y = r3 * vtx[i3].p.y;
      elm[i][j].p3.z = r3 * vtx[i3].p.z;

      /* Set coordinates for vertex 4 of the element using current shell radius */
      elm[i][j].p4.x = r4 * vtx[i1].p.x;
      elm[i][j].p4.y = r4 * vtx[i1].p.y;
      elm[i][j].p4.z = r4 * vtx[i1].p.z;

      /* Set coordinates for vertex 5 of the element using current shell radius */
      elm[i][j].p5.x = r5 * vtx[i2].p.x;
      elm[i][j].p5.y = r5 * vtx[i2].p.y;
      elm[i][j].p5.z = r5 * vtx[i2].p.z;

      /* Set coordinates for vertex 6 of the element using current shell radius */
      elm[i][j].p6.x = r6 * vtx[i3].p.x;
      elm[i][j].p6.y = r6 * vtx[i3].p.y;
      elm[i][j].p6.z = r6 * vtx[i3].p.z;

      /* Set density for all vertices using previous shell value converted to kg/km^3 */
      elm[i][j].rho1 = rho[i1][i - 1] * 1E12;
      elm[i][j].rho2 = rho[i2][i - 1] * 1E12;
      elm[i][j].rho3 = rho[i3][i - 1] * 1E12;
      elm[i][j].rho4 = rho[i1][i]     * 1E12;
      elm[i][j].rho5 = rho[i2][i]     * 1E12;
      elm[i][j].rho6 = rho[i3][i]     * 1E12;
    }
  }
}

static void decomposeElements (int nl, int nf_l,
                               struct Element elm[nl][nf_l],
                               struct Prism prism[nl][nf_l][N_TETRAHEDRA])
{
  /* Decomposes prismatic elements into tetrahedra with computed densities and
     centers of mass. */
  for (int i = 0; i < nl; i++)

    for (int j = 0; j < nf_l; j++)
    {
      /* Extract density values from the current element */
      double rho1 = elm[i][j].rho1;
      double rho2 = elm[i][j].rho2;
      double rho3 = elm[i][j].rho3;
      double rho4 = elm[i][j].rho4;
      double rho5 = elm[i][j].rho5;
      double rho6 = elm[i][j].rho6;

      /* Compute the average density of the element */
      double rhom = (rho1 + rho2 + rho3 +
                     rho4 + rho5 + rho6) / 6;

      /* Calculate average densities between adjacent vertices */
      double rho12 = 0.5 * (rho1 + rho2);
      double rho13 = 0.5 * (rho1 + rho3);
      double rho23 = 0.5 * (rho2 + rho3);
      double rho45 = 0.5 * (rho4 + rho5);
      double rho46 = 0.5 * (rho4 + rho6);
      double rho56 = 0.5 * (rho5 + rho6);
      double rho14 = 0.5 * (rho1 + rho4);
      double rho25 = 0.5 * (rho2 + rho5);
      double rho36 = 0.5 * (rho3 + rho6);

      /* Calculate average densities for quadrilateral faces */
      double rho1245 = 0.25 * (rho1 + rho2 + rho4 + rho5);
      double rho2356 = 0.25 * (rho2 + rho3 + rho5 + rho6);
      double rho1346 = 0.25 * (rho1 + rho3 + rho4 + rho6);

      struct Triplet p1, p2, p3, p4, p5, p6,
                   pm, p12, p13, p23, p45,
                   p46, p56, p14, p25, p36,
                   p1245, p2356, p1346;

      /* Assign coordinates for vertex 1 */
      p1.x = elm[i][j].p1.x;
      p1.y = elm[i][j].p1.y;
      p1.z = elm[i][j].p1.z;

      /* Assign coordinates for vertex 2 */
      p2.x = elm[i][j].p2.x;
      p2.y = elm[i][j].p2.y;
      p2.z = elm[i][j].p2.z;

      /* Assign coordinates for vertex 3 */
      p3.x = elm[i][j].p3.x;
      p3.y = elm[i][j].p3.y;
      p3.z = elm[i][j].p3.z;

      /* Assign coordinates for vertex 4 */
      p4.x = elm[i][j].p4.x;
      p4.y = elm[i][j].p4.y;
      p4.z = elm[i][j].p4.z;

      /* Assign coordinates for vertex 5 */
      p5.x = elm[i][j].p5.x;
      p5.y = elm[i][j].p5.y;
      p5.z = elm[i][j].p5.z;

      /* Assign coordinates for vertex 6 */
      p6.x = elm[i][j].p6.x;
      p6.y = elm[i][j].p6.y;
      p6.z = elm[i][j].p6.z;

      /* Compute the centroid of the element */
      pm.x = (p1.x + p2.x + p3.x +
              p4.x + p5.x + p6.x) / 6;
      pm.y = (p1.y + p2.y + p3.y +
              p4.y + p5.y + p6.y) / 6;
      pm.z = (p1.z + p2.z + p3.z +
              p4.z + p5.z + p6.z) / 6;

      /* Compute midpoint coordinates between vertex 1 and 2 */
      p12.x = 0.5 * (p1.x + p2.x);
      p12.y = 0.5 * (p1.y + p2.y);
      p12.z = 0.5 * (p1.z + p2.z);

      /* Compute midpoint coordinates between vertex 1 and 3 */
      p13.x = 0.5 * (p1.x + p3.x);
      p13.y = 0.5 * (p1.y + p3.y);
      p13.z = 0.5 * (p1.z + p3.z);

      /* Compute midpoint coordinates between vertex 2 and 3 */
      p23.x = 0.5 * (p2.x + p3.x);
      p23.y = 0.5 * (p2.y + p3.y);
      p23.z = 0.5 * (p2.z + p3.z);

      /* Compute midpoint coordinates between vertex 4 and 5 */
      p45.x = 0.5 * (p4.x + p5.x);
      p45.y = 0.5 * (p4.y + p5.y);
      p45.z = 0.5 * (p4.z + p5.z);

      /* Compute midpoint coordinates between vertex 4 and 6 */
      p46.x = 0.5 * (p4.x + p6.x);
      p46.y = 0.5 * (p4.y + p6.y);
      p46.z = 0.5 * (p4.z + p6.z);

      /* Compute midpoint coordinates between vertex 5 and 6 */
      p56.x = 0.5 * (p5.x + p6.x);
      p56.y = 0.5 * (p5.y + p6.y);
      p56.z = 0.5 * (p5.z + p6.z);

      /* Compute midpoint coordinates between vertex 1 and 4 */
      p14.x = 0.5 * (p1.x + p4.x);
      p14.y = 0.5 * (p1.y + p4.y);
      p14.z = 0.5 * (p1.z + p4.z);

      /* Compute midpoint coordinates between vertex 2 and 5 */
      p25.x = 0.5 * (p2.x + p5.x);
      p25.y = 0.5 * (p2.y + p5.y);
      p25.z = 0.5 * (p2.z + p5.z);

      /* Compute midpoint coordinates between vertex 3 and 6 */
      p36.x = 0.5 * (p3.x + p6.x);
      p36.y = 0.5 * (p3.y + p6.y);
      p36.z = 0.5 * (p3.z + p6.z);

      /* Compute centroid coordinates for quadrilateral face 1-2-4-5 */
      p1245.x = 0.25 * (p1.x + p2.x + p4.x + p5.x);
      p1245.y = 0.25 * (p1.y + p2.y + p4.y + p5.y);
      p1245.z = 0.25 * (p1.z + p2.z + p4.z + p5.z);

      /* Compute centroid coordinates for quadrilateral face 2-3-5-6 */
      p2356.x = 0.25 * (p2.x + p3.x + p5.x + p6.x);
      p2356.y = 0.25 * (p2.y + p3.y + p5.y + p6.y);
      p2356.z = 0.25 * (p2.z + p3.z + p5.z + p6.z);

      /* Compute centroid coordinates for quadrilateral face 1-3-4-6 */
      p1346.x = 0.25 * (p1.x + p3.x + p4.x + p6.x);
      p1346.y = 0.25 * (p1.y + p3.y + p4.y + p6.y);
      p1346.z = 0.25 * (p1.z + p3.z + p4.z + p6.z);

      struct Tetrahedron th1, th2, th3, th4, th5, th6,
                         th7, th8, th9, th10, th11, th12,
                         th13, th14, th15, th16, th17,
                         th18, th19, th20, th21, th22,
                         th23, th24, th25, th26, th27,
                         th28, th29, th30, th31, th32;

      /* Define tetrahedrons 1 to 4 with centroid and vertices */
      th1.p1  = pm; th1.p2 = p1;  th1.p3 = p12; th1.p4 = p13;
      th2.p1  = pm; th2.p2 = p2;  th2.p3 = p12; th2.p4 = p23;
      th3.p1  = pm; th3.p2 = p3;  th3.p3 = p13; th3.p4 = p23;
      th4.p1  = pm; th4.p2 = p12; th4.p3 = p13; th4.p4 = p23;

      /* Define tetrahedrons 5 to 8 with centroid and vertices */
      th5.p1  = pm; th5.p2 = p4;  th5.p3 = p45; th5.p4 = p46;
      th6.p1  = pm; th6.p2 = p5;  th6.p3 = p45; th6.p4 = p56;
      th7.p1  = pm; th7.p2 = p6;  th7.p3 = p46; th7.p4 = p56;
      th8.p1  = pm; th8.p2 = p45; th8.p3 = p46; th8.p4 = p56;

      /* Define tetrahedrons 9 to 25 with centroid and vertices */
      th9.p1  = pm; th9.p2  = p1;    th9.p3  = p12;   th9.p4  = p14;
      th10.p1 = pm; th10.p2 = p1245; th10.p3 = p12;   th10.p4 = p14;
      th11.p1 = pm; th11.p2 = p12;   th11.p3 = p2;    th11.p4 = p25;
      th12.p1 = pm; th12.p2 = p12;   th12.p3 = p1245; th12.p4 = p25;
      th13.p1 = pm; th13.p2 = p14;   th13.p3 = p4;    th13.p4 = p45;
      th14.p1 = pm; th14.p2 = p14;   th14.p3 = p1245; th14.p4 = p45;
      th15.p1 = pm; th15.p2 = p25;   th15.p3 = p5;    th15.p4 = p45;
      th16.p1 = pm; th16.p2 = p25;   th16.p3 = p1245; th16.p4 = p45;

      /* Define tetrahedrons 17 to 24 with centroid and vertices */
      th17.p1 = pm; th17.p2 = p2;    th17.p3 = p23;   th17.p4 = p25;
      th18.p1 = pm; th18.p2 = p2356; th18.p3 = p23;   th18.p4 = p25;
      th19.p1 = pm; th19.p2 = p23;   th19.p3 = p3;    th19.p4 = p36;
      th20.p1 = pm; th20.p2 = p23;   th20.p3 = p2356; th20.p4 = p36;
      th21.p1 = pm; th21.p2 = p5;    th21.p3 = p25;   th21.p4 = p56;
      th22.p1 = pm; th22.p2 = p2356; th22.p3 = p25;   th22.p4 = p56;
      th23.p1 = pm; th23.p2 = p36;   th23.p3 = p56;   th23.p4 = p6;
      th24.p1 = pm; th24.p2 = p36;   th24.p3 = p56;   th24.p4 = p2356;

      /* Define tetrahedrons 25 to 32 with centroid and vertices */
      th25.p1 = pm; th25.p2 = p1;    th25.p3 = p13;   th25.p4 = p14;
      th26.p1 = pm; th26.p2 = p1346; th26.p3 = p13;   th26.p4 = p14;
      th27.p1 = pm; th27.p2 = p13;   th27.p3 = p3;    th27.p4 = p36;
      th28.p1 = pm; th28.p2 = p13;   th28.p3 = p1346; th28.p4 = p36;
      th29.p1 = pm; th29.p2 = p14;   th29.p3 = p4;    th29.p4 = p46;
      th30.p1 = pm; th30.p2 = p14;   th30.p3 = p1346; th30.p4 = p46;
      th31.p1 = pm; th31.p2 = p36;   th31.p3 = p46;   th31.p4 = p6;
      th32.p1 = pm; th32.p2 = p36;   th32.p3 = p46;   th32.p4 = p1346;

      /* Compute sum of densities for tetrahedrons 1 to 4 */
      double sum_th1  = rhom + rho1  + rho12 + rho13;
      double sum_th2  = rhom + rho2  + rho12 + rho23;
      double sum_th3  = rhom + rho3  + rho13 + rho23;
      double sum_th4  = rhom + rho12 + rho13 + rho23;

      /* Compute sum of densities for tetrahedrons 5 to 8 */
      double sum_th5  = rhom + rho4  + rho45 + rho46;
      double sum_th6  = rhom + rho5  + rho45 + rho56;
      double sum_th7  = rhom + rho6  + rho46 + rho56;
      double sum_th8  = rhom + rho45 + rho46 + rho56;

      /* Compute sum of densities for tetrahedrons 9 to 16 */
      double sum_th9  = rhom + rho1    + rho12   + rho14;
      double sum_th10 = rhom + rho1245 + rho12   + rho14;
      double sum_th11 = rhom + rho12   + rho2    + rho25;
      double sum_th12 = rhom + rho12   + rho1245 + rho25;
      double sum_th13 = rhom + rho14   + rho4    + rho45;
      double sum_th14 = rhom + rho14   + rho1245 + rho45;
      double sum_th15 = rhom + rho25   + rho5    + rho45;
      double sum_th16 = rhom + rho25   + rho1245 + rho45;

      /* Compute sum of densities for tetrahedrons 17 to 24 */
      double sum_th17 = rhom + rho2    + rho23   + rho25;
      double sum_th18 = rhom + rho2356 + rho23   + rho25;
      double sum_th19 = rhom + rho23   + rho3    + rho36;
      double sum_th20 = rhom + rho23   + rho2356 + rho36;
      double sum_th21 = rhom + rho5    + rho25   + rho56;
      double sum_th22 = rhom + rho2356 + rho25   + rho56;
      double sum_th23 = rhom + rho36   + rho56   + rho6;
      double sum_th24 = rhom + rho36   + rho56   + rho2356;

      /* Compute sum of densities for tetrahedrons 25 to 32 */
      double sum_th25 = rhom + rho1    + rho13   + rho14;
      double sum_th26 = rhom + rho1346 + rho13   + rho14;
      double sum_th27 = rhom + rho13   + rho3    + rho36;
      double sum_th28 = rhom + rho13   + rho1346 + rho36;
      double sum_th29 = rhom + rho14   + rho4    + rho46;
      double sum_th30 = rhom + rho14   + rho1346 + rho46;
      double sum_th31 = rhom + rho36   + rho46   + rho6;
      double sum_th32 = rhom + rho36   + rho46   + rho1346;

      /* Compute average density for tetrahedrons 1 to 4 */
      double rho_th1  = 0.25 * sum_th1;
      double rho_th2  = 0.25 * sum_th2;
      double rho_th3  = 0.25 * sum_th3;
      double rho_th4  = 0.25 * sum_th4;

      /* Compute average density for tetrahedrons 5 to 8 */
      double rho_th5  = 0.25 * sum_th5;
      double rho_th6  = 0.25 * sum_th6;
      double rho_th7  = 0.25 * sum_th7;
      double rho_th8  = 0.25 * sum_th8;

      /* Compute average density for tetrahedrons 9 to 16 */
      double rho_th9  = 0.25 * sum_th9;
      double rho_th10 = 0.25 * sum_th10;
      double rho_th11 = 0.25 * sum_th11;
      double rho_th12 = 0.25 * sum_th12;
      double rho_th13 = 0.25 * sum_th13;
      double rho_th14 = 0.25 * sum_th14;
      double rho_th15 = 0.25 * sum_th15;
      double rho_th16 = 0.25 * sum_th16;

      /* Compute average density for tetrahedrons 17 to 24 */
      double rho_th17 = 0.25 * sum_th17;
      double rho_th18 = 0.25 * sum_th18;
      double rho_th19 = 0.25 * sum_th19;
      double rho_th20 = 0.25 * sum_th20;
      double rho_th21 = 0.25 * sum_th21;
      double rho_th22 = 0.25 * sum_th22;
      double rho_th23 = 0.25 * sum_th23;
      double rho_th24 = 0.25 * sum_th24;

      /* Compute average density for tetrahedrons 25 to 32 */
      double rho_th25 = 0.25 * sum_th25;
      double rho_th26 = 0.25 * sum_th26;
      double rho_th27 = 0.25 * sum_th27;
      double rho_th28 = 0.25 * sum_th28;
      double rho_th29 = 0.25 * sum_th29;
      double rho_th30 = 0.25 * sum_th30;
      double rho_th31 = 0.25 * sum_th31;
      double rho_th32 = 0.25 * sum_th32;

      struct Triplet cm1, cm2, cm3, cm4, cm5, cm6,
                   cm7, cm8, cm9, cm10, cm11, cm12,
                   cm13, cm14, cm15, cm16, cm17,
                   cm18, cm19, cm20, cm21, cm22,
                   cm23, cm24, cm25, cm26, cm27,
                   cm28, cm29, cm30, cm31, cm32;

      /* Compute center of mass for tetrahedron 1 */
      cm1.x  = (rhom  * pm.x  + rho1  * p1.x   +
                rho12 * p12.x + rho13 * p13.x) / (sum_th1 + EPSILON);
      cm1.y  = (rhom  * pm.y  + rho1  * p1.y   +
                rho12 * p12.y + rho13 * p13.y) / (sum_th1 + EPSILON);
      cm1.z  = (rhom  * pm.z  + rho1  * p1.z   +
                rho12 * p12.z + rho13 * p13.z) / (sum_th1 + EPSILON);

      /* Compute center of mass for tetrahedron 2 */
      cm2.x  = (rhom  * pm.x  + rho2  * p2.x   +
                rho12 * p12.x + rho23 * p23.x) / (sum_th2 + EPSILON);
      cm2.y  = (rhom  * pm.y  + rho2  * p2.y   +
                rho12 * p12.y + rho23 * p23.y) / (sum_th2 + EPSILON);
      cm2.z  = (rhom  * pm.z  + rho2  * p2.z   +
                rho12 * p12.z + rho23 * p23.z) / (sum_th2 + EPSILON);

      /* Compute center of mass for tetrahedron 3 */
      cm3.x  = (rhom  * pm.x  + rho3  * p3.x   +
                rho13 * p13.x + rho23 * p23.x) / (sum_th3 + EPSILON);
      cm3.y  = (rhom  * pm.y  + rho3  * p3.y   +
                rho13 * p13.y + rho23 * p23.y) / (sum_th3 + EPSILON);
      cm3.z  = (rhom  * pm.z  + rho3  * p3.z   +
                rho13 * p13.z + rho23 * p23.z) / (sum_th3 + EPSILON);

      /* Compute center of mass for tetrahedron 4 */
      cm4.x  = (rhom  * pm.x  + rho12 * p12.x  +
                rho13 * p13.x + rho23 * p23.x) / (sum_th4 + EPSILON);
      cm4.y  = (rhom  * pm.y  + rho12 * p12.y  +
                rho13 * p13.y + rho23 * p23.y) / (sum_th4 + EPSILON);
      cm4.z  = (rhom  * pm.z  + rho12 * p12.z  +
                rho13 * p13.z + rho23 * p23.z) / (sum_th4 + EPSILON);

      /* Compute center of mass for tetrahedron 5 */
      cm5.x  = (rhom  * pm.x  + rho4  * p4.x   +
                rho45 * p45.x + rho46 * p46.x) / (sum_th5 + EPSILON);
      cm5.y  = (rhom  * pm.y  + rho4  * p4.y   +
                rho45 * p45.y + rho46 * p46.y) / (sum_th5 + EPSILON);
      cm5.z  = (rhom  * pm.z  + rho4  * p4.z   +
                rho45 * p45.z + rho46 * p46.z) / (sum_th5 + EPSILON);

      /* Compute center of mass for tetrahedron 6 */
      cm6.x  = (rhom  * pm.x  + rho5  * p5.x   +
                rho45 * p45.x + rho56 * p56.x) / (sum_th6 + EPSILON);
      cm6.y  = (rhom  * pm.y  + rho5  * p5.y   +
                rho45 * p45.y + rho56 * p56.y) / (sum_th6 + EPSILON);
      cm6.z  = (rhom  * pm.z  + rho5  * p5.z   +
                rho45 * p45.z + rho56 * p56.z) / (sum_th6 + EPSILON);

      /* Compute center of mass for tetrahedron 7 */
      cm7.x  = (rhom  * pm.x  + rho6  * p6.x   +
                rho46 * p46.x + rho56 * p56.x) / (sum_th7 + EPSILON);
      cm7.y  = (rhom  * pm.y  + rho6  * p6.y   +
                rho46 * p46.y + rho56 * p56.y) / (sum_th7 + EPSILON);
      cm7.z  = (rhom  * pm.z  + rho6  * p6.z   +
                rho46 * p46.z + rho56 * p56.z) / (sum_th7 + EPSILON);

      /* Compute center of mass for tetrahedron 8 */
      cm8.x  = (rhom  * pm.x  + rho45 * p45.x  +
                rho46 * p46.x + rho56 * p56.x) / (sum_th8 + EPSILON);
      cm8.y  = (rhom  * pm.y  + rho45 * p45.y  +
                rho46 * p46.y + rho56 * p56.y) / (sum_th8 + EPSILON);
      cm8.z  = (rhom  * pm.z  + rho45 * p45.z  +
                rho46 * p46.z + rho56 * p56.z) / (sum_th8 + EPSILON);

      /* Compute center of mass for tetrahedron 9 */
      cm9.x  = (rhom    * pm.x    + rho1    * p1.x    +
                rho12   * p12.x   + rho14   * p14.x)  / (sum_th9 + EPSILON);
      cm9.y  = (rhom    * pm.y    + rho1    * p1.y    +
                rho12   * p12.y   + rho14   * p14.y)  / (sum_th9 + EPSILON);
      cm9.z  = (rhom    * pm.z    + rho1    * p1.z    +
                rho12   * p12.z   + rho14   * p14.z)  / (sum_th9 + EPSILON);

      /* Compute center of mass for tetrahedron 10 */
      cm10.x = (rhom    * pm.x    + rho1245 * p1245.x +
                rho12   * p12.x   + rho14   * p14.x)  / (sum_th10 + EPSILON);
      cm10.y = (rhom    * pm.y    + rho1245 * p1245.y +
                rho12   * p12.y   + rho14   * p14.y)  / (sum_th10 + EPSILON);
      cm10.z = (rhom    * pm.z    + rho1245 * p1245.z +
                rho12   * p12.z   + rho14   * p14.z)  / (sum_th10 + EPSILON);

      /* Compute center of mass for tetrahedron 11 */
      cm11.x = (rhom    * pm.x    + rho12   * p12.x   +
                rho2    * p2.x    + rho25   * p25.x)  / (sum_th11 + EPSILON);
      cm11.y = (rhom    * pm.y    + rho12   * p12.y   +
                rho2    * p2.y    + rho25   * p25.y)  / (sum_th11 + EPSILON);
      cm11.z = (rhom    * pm.z    + rho12   * p12.z   +
                rho2    * p2.z    + rho25   * p25.z)  / (sum_th11 + EPSILON);

      /* Compute center of mass for tetrahedron 12 */
      cm12.x = (rhom    * pm.x    + rho12   * p12.x   +
                rho1245 * p1245.x + rho25   * p25.x)  / (sum_th12 + EPSILON);
      cm12.y = (rhom    * pm.y    + rho12   * p12.y   +
                rho1245 * p1245.y + rho25   * p25.y)  / (sum_th12 + EPSILON);
      cm12.z = (rhom    * pm.z    + rho12   * p12.z   +
                rho1245 * p1245.z + rho25   * p25.z)  / (sum_th12 + EPSILON);

      /* Compute center of mass for tetrahedron 13 */
      cm13.x = (rhom    * pm.x    + rho14   * p14.x   +
                rho4    * p4.x    + rho45   * p45.x)  / (sum_th13 + EPSILON);
      cm13.y = (rhom    * pm.y    + rho14   * p14.y   +
                rho4    * p4.y    + rho45   * p45.y)  / (sum_th13 + EPSILON);
      cm13.z = (rhom    * pm.z    + rho14   * p14.z   +
                rho4    * p4.z    + rho45   * p45.z)  / (sum_th13 + EPSILON);

      /* Compute center of mass for tetrahedron 14 */
      cm14.x = (rhom    * pm.x    + rho14   * p14.x   +
                rho1245 * p1245.x + rho45   * p45.x)  / (sum_th14 + EPSILON);
      cm14.y = (rhom    * pm.y    + rho14   * p14.y   +
                rho1245 * p1245.y + rho45   * p45.y)  / (sum_th14 + EPSILON);
      cm14.z = (rhom    * pm.z    + rho14   * p14.z   +
                rho1245 * p1245.z + rho45   * p45.z)  / (sum_th14 + EPSILON);

      /* Compute center of mass for tetrahedron 15 */
      cm15.x = (rhom    * pm.x    + rho25   * p25.x   +
                rho5    * p5.x    + rho45   * p45.x)  / (sum_th15 + EPSILON);
      cm15.y = (rhom    * pm.y    + rho25   * p25.y   +
                rho5    * p5.y    + rho45   * p45.y)  / (sum_th15 + EPSILON);
      cm15.z = (rhom    * pm.z    + rho25   * p25.z   +
                rho5    * p5.z    + rho45   * p45.z)  / (sum_th15 + EPSILON);

      /* Compute center of mass for tetrahedron 16 */
      cm16.x = (rhom    * pm.x    + rho25   * p25.x   +
                rho1245 * p1245.x + rho45   * p45.x)  / (sum_th16 + EPSILON);
      cm16.y = (rhom    * pm.y    + rho25   * p25.y   +
                rho1245 * p1245.y + rho45   * p45.y)  / (sum_th16 + EPSILON);
      cm16.z = (rhom    * pm.z    + rho25   * p25.z   +
                rho1245 * p1245.z + rho45   * p45.z)  / (sum_th16 + EPSILON);

      /* Compute center of mass for tetrahedron 17 */
      cm17.x = (rhom    * pm.x    + rho2    * p2.x     +
                rho23   * p23.x   + rho25   * p25.x)   / (sum_th17 + EPSILON);
      cm17.y = (rhom    * pm.y    + rho2    * p2.y     +
                rho23   * p23.y   + rho25   * p25.y)   / (sum_th17 + EPSILON);
      cm17.z = (rhom    * pm.z    + rho2    * p2.z     +
                rho23   * p23.z   + rho25   * p25.z)   / (sum_th17 + EPSILON);

      /* Compute center of mass for tetrahedron 18 */
      cm18.x = (rhom    * pm.x    + rho2356 * p2356.x  +
                rho23   * p23.x   + rho25   * p25.x)   / (sum_th18 + EPSILON);
      cm18.y = (rhom    * pm.y    + rho2356 * p2356.y  +
                rho23   * p23.y   + rho25   * p25.y)   / (sum_th18 + EPSILON);
      cm18.z = (rhom    * pm.z    + rho2356 * p2356.z  +
                rho23   * p23.z   + rho25   * p25.z)   / (sum_th18 + EPSILON);

      /* Compute center of mass for tetrahedron 19 */
      cm19.x = (rhom    * pm.x    + rho23   * p23.x    +
                rho3    * p3.x    + rho36   * p36.x)   / (sum_th19 + EPSILON);
      cm19.y = (rhom    * pm.y    + rho23   * p23.y    +
                rho3    * p3.y    + rho36   * p36.y)   / (sum_th19 + EPSILON);
      cm19.z = (rhom    * pm.z    + rho23   * p23.z    +
                rho3    * p3.z    + rho36   * p36.z)   / (sum_th19 + EPSILON);

      /* Compute center of mass for tetrahedron 20 */
      cm20.x = (rhom    * pm.x    + rho23   * p23.x    +
                rho2356 * p2356.x + rho36   * p36.x)   / (sum_th20 + EPSILON);
      cm20.y = (rhom    * pm.y    + rho23   * p23.y    +
                rho2356 * p2356.y + rho36   * p36.y)   / (sum_th20 + EPSILON);
      cm20.z = (rhom    * pm.z    + rho23   * p23.z    +
                rho2356 * p2356.z + rho36   * p36.z)   / (sum_th20 + EPSILON);

      /* Compute center of mass for tetrahedron 21 */
      cm21.x = (rhom    * pm.x    + rho5    * p5.x     +
                rho25   * p25.x   + rho56   * p56.x)   / (sum_th21 + EPSILON);
      cm21.y = (rhom    * pm.y    + rho5    * p5.y     +
                rho25   * p25.y   + rho56   * p56.y)   / (sum_th21 + EPSILON);
      cm21.z = (rhom    * pm.z    + rho5    * p5.z     +
                rho25   * p25.z   + rho56   * p56.z)   / (sum_th21 + EPSILON);

      /* Compute center of mass for tetrahedron 22 */
      cm22.x = (rhom    * pm.x    + rho2356 * p2356.x  +
                rho25   * p25.x   + rho56   * p56.x)   / (sum_th22 + EPSILON);
      cm22.y = (rhom    * pm.y    + rho2356 * p2356.y  +
                rho25   * p25.y   + rho56   * p56.y)   / (sum_th22 + EPSILON);
      cm22.z = (rhom    * pm.z    + rho2356 * p2356.z  +
                rho25   * p25.z   + rho56   * p56.z)   / (sum_th22 + EPSILON);

      /* Compute center of mass for tetrahedron 23 */
      cm23.x = (rhom    * pm.x    + rho36   * p36.x    +
                rho56   * p56.x   + rho6    * p6.x)    / (sum_th23 + EPSILON);
      cm23.y = (rhom    * pm.y    + rho36   * p36.y    +
                rho56   * p56.y   + rho6    * p6.y)    / (sum_th23 + EPSILON);
      cm23.z = (rhom    * pm.z    + rho36   * p36.z    +
                rho56   * p56.z   + rho6    * p6.z)    / (sum_th23 + EPSILON);

      /* Compute center of mass for tetrahedron 24 */
      cm24.x = (rhom    * pm.x    + rho36   * p36.x    +
                rho56   * p56.x   + rho2356 * p2356.x) / (sum_th24 + EPSILON);
      cm24.y = (rhom    * pm.y    + rho36   * p36.y    +
                rho56   * p56.y   + rho2356 * p2356.y) / (sum_th24 + EPSILON);
      cm24.z = (rhom    * pm.z    + rho36   * p36.z    +
                rho56   * p56.z   + rho2356 * p2356.z) / (sum_th24 + EPSILON);

      /* Compute center of mass for tetrahedron 25 */
      cm25.x = (rhom    * pm.x    + rho1    * p1.x     +
                rho13   * p13.x   + rho14   * p14.x)   / (sum_th25 + EPSILON);
      cm25.y = (rhom    * pm.y    + rho1    * p1.y     +
                rho13   * p13.y   + rho14   * p14.y)   / (sum_th25 + EPSILON);
      cm25.z = (rhom    * pm.z    + rho1    * p1.z     +
                rho13   * p13.z   + rho14   * p14.z)   / (sum_th25 + EPSILON);

      /* Compute center of mass for tetrahedron 26 */
      cm26.x = (rhom    * pm.x    + rho1346 * p1346.x  +
                rho13   * p13.x   + rho14   * p14.x)   / (sum_th26 + EPSILON);
      cm26.y = (rhom    * pm.y    + rho1346 * p1346.y  +
                rho13   * p13.y   + rho14   * p14.y)   / (sum_th26 + EPSILON);
      cm26.z = (rhom    * pm.z    + rho1346 * p1346.z  +
                rho13   * p13.z   + rho14   * p14.z)   / (sum_th26 + EPSILON);

      /* Compute center of mass for tetrahedron 27 */
      cm27.x = (rhom    * pm.x    + rho13   * p13.x    +
                rho3    * p3.x    + rho36   * p36.x)   / (sum_th27 + EPSILON);
      cm27.y = (rhom    * pm.y    + rho13   * p13.y    +
                rho3    * p3.y    + rho36   * p36.y)   / (sum_th27 + EPSILON);
      cm27.z = (rhom    * pm.z    + rho13   * p13.z    +
                rho3    * p3.z    + rho36   * p36.z)   / (sum_th27 + EPSILON);

      /* Compute center of mass for tetrahedron 28 */
      cm28.x = (rhom    * pm.x    + rho13   * p13.x    +
                rho1346 * p1346.x + rho36   * p36.x)   / (sum_th28 + EPSILON);
      cm28.y = (rhom    * pm.y    + rho13   * p13.y    +
                rho1346 * p1346.y + rho36   * p36.y)   / (sum_th28 + EPSILON);
      cm28.z = (rhom    * pm.z    + rho13   * p13.z    +
                rho1346 * p1346.z + rho36   * p36.z)   / (sum_th28 + EPSILON);

      /* Compute center of mass for tetrahedron 29 */
      cm29.x = (rhom    * pm.x    + rho14   * p14.x    +
                rho4    * p4.x    + rho46   * p46.x)   / (sum_th29 + EPSILON);
      cm29.y = (rhom    * pm.y    + rho14   * p14.y    +
                rho4    * p4.y    + rho46   * p46.y)   / (sum_th29 + EPSILON);
      cm29.z = (rhom    * pm.z    + rho14   * p14.z    +
                rho4    * p4.z    + rho46   * p46.z)   / (sum_th29 + EPSILON);

      /* Compute center of mass for tetrahedron 30 */
      cm30.x = (rhom    * pm.x    + rho14   * p14.x    +
                rho1346 * p1346.x + rho46   * p46.x)   / (sum_th30 + EPSILON);
      cm30.y = (rhom    * pm.y    + rho14   * p14.y    +
                rho1346 * p1346.y + rho46   * p46.y)   / (sum_th30 + EPSILON);
      cm30.z = (rhom    * pm.z    + rho14   * p14.z    +
                rho1346 * p1346.z + rho46   * p46.z)   / (sum_th30 + EPSILON);

      /* Compute center of mass for tetrahedron 31 */
      cm31.x = (rhom    * pm.x    + rho36   * p36.x    +
                rho46   * p46.x   + rho6    * p6.x)    / (sum_th31 + EPSILON);
      cm31.y = (rhom    * pm.y    + rho36   * p36.y    +
                rho46   * p46.y   + rho6    * p6.y)    / (sum_th31 + EPSILON);
      cm31.z = (rhom    * pm.z    + rho36   * p36.z    +
                rho46   * p46.z   + rho6    * p6.z)    / (sum_th31 + EPSILON);

      /* Compute center of mass for tetrahedron 32 */
      cm32.x = (rhom    * pm.x    + rho36   * p36.x    +
                rho46   * p46.x   + rho1346 * p1346.x) / (sum_th32 + EPSILON);
      cm32.y = (rhom    * pm.y    + rho36   * p36.y    +
                rho46   * p46.y   + rho1346 * p1346.y) / (sum_th32 + EPSILON);
      cm32.z = (rhom    * pm.z    + rho36   * p36.z    +
                rho46   * p46.z   + rho1346 * p1346.z) / (sum_th32 + EPSILON);

      /* Compute volume for tetrahedrons 1 to 4 */
      double dv1  = volumeTetrahedron (&th1);
      double dv2  = volumeTetrahedron (&th2);
      double dv3  = volumeTetrahedron (&th3);
      double dv4  = volumeTetrahedron (&th4);

      /* Compute volume for tetrahedrons 5 to 8 */
      double dv5  = volumeTetrahedron (&th5);
      double dv6  = volumeTetrahedron (&th6);
      double dv7  = volumeTetrahedron (&th7);
      double dv8  = volumeTetrahedron (&th8);

      /* Compute volume for tetrahedrons 9 to 16 */
      double dv9  = volumeTetrahedron (&th9);
      double dv10 = volumeTetrahedron (&th10);
      double dv11 = volumeTetrahedron (&th11);
      double dv12 = volumeTetrahedron (&th12);
      double dv13 = volumeTetrahedron (&th13);
      double dv14 = volumeTetrahedron (&th14);
      double dv15 = volumeTetrahedron (&th15);
      double dv16 = volumeTetrahedron (&th16);

      /* Compute volume for tetrahedrons 17 to 24 */
      double dv17 = volumeTetrahedron (&th17);
      double dv18 = volumeTetrahedron (&th18);
      double dv19 = volumeTetrahedron (&th19);
      double dv20 = volumeTetrahedron (&th20);
      double dv21 = volumeTetrahedron (&th21);
      double dv22 = volumeTetrahedron (&th22);
      double dv23 = volumeTetrahedron (&th23);
      double dv24 = volumeTetrahedron (&th24);

      /* Compute volume for tetrahedrons 25 to 32 */
      double dv25 = volumeTetrahedron (&th25);
      double dv26 = volumeTetrahedron (&th26);
      double dv27 = volumeTetrahedron (&th27);
      double dv28 = volumeTetrahedron (&th28);
      double dv29 = volumeTetrahedron (&th29);
      double dv30 = volumeTetrahedron (&th30);
      double dv31 = volumeTetrahedron (&th31);
      double dv32 = volumeTetrahedron (&th32);

      /* Set skip flag for tetrahedrons 1 to 4 based on volume threshold */
      prism[i][j][0].skip  = dv1  < SKIP_VOLUME ? true : false;
      prism[i][j][1].skip  = dv2  < SKIP_VOLUME ? true : false;
      prism[i][j][2].skip  = dv3  < SKIP_VOLUME ? true : false;
      prism[i][j][3].skip  = dv4  < SKIP_VOLUME ? true : false;

      /* Set skip flag for tetrahedrons 5 to 8 based on volume threshold */
      prism[i][j][4].skip  = dv5  < SKIP_VOLUME ? true : false;
      prism[i][j][5].skip  = dv6  < SKIP_VOLUME ? true : false;
      prism[i][j][6].skip  = dv7  < SKIP_VOLUME ? true : false;
      prism[i][j][7].skip  = dv8  < SKIP_VOLUME ? true : false;

      /* Set skip flag for tetrahedrons 9 to 16 based on volume threshold */
      prism[i][j][8].skip  = dv9  < SKIP_VOLUME ? true : false;
      prism[i][j][9].skip  = dv10 < SKIP_VOLUME ? true : false;
      prism[i][j][10].skip = dv11 < SKIP_VOLUME ? true : false;
      prism[i][j][11].skip = dv12 < SKIP_VOLUME ? true : false;
      prism[i][j][12].skip = dv13 < SKIP_VOLUME ? true : false;
      prism[i][j][13].skip = dv14 < SKIP_VOLUME ? true : false;
      prism[i][j][14].skip = dv15 < SKIP_VOLUME ? true : false;
      prism[i][j][15].skip = dv16 < SKIP_VOLUME ? true : false;

      /* Set skip flag for tetrahedrons 17 to 24 based on volume threshold */
      prism[i][j][16].skip = dv17 < SKIP_VOLUME ? true : false;
      prism[i][j][17].skip = dv18 < SKIP_VOLUME ? true : false;
      prism[i][j][18].skip = dv19 < SKIP_VOLUME ? true : false;
      prism[i][j][19].skip = dv20 < SKIP_VOLUME ? true : false;
      prism[i][j][20].skip = dv21 < SKIP_VOLUME ? true : false;
      prism[i][j][21].skip = dv22 < SKIP_VOLUME ? true : false;
      prism[i][j][22].skip = dv23 < SKIP_VOLUME ? true : false;
      prism[i][j][23].skip = dv24 < SKIP_VOLUME ? true : false;

      /* Set skip flag for tetrahedrons 25 to 32 based on volume threshold */
      prism[i][j][24].skip = dv25 < SKIP_VOLUME ? true : false;
      prism[i][j][25].skip = dv26 < SKIP_VOLUME ? true : false;
      prism[i][j][26].skip = dv27 < SKIP_VOLUME ? true : false;
      prism[i][j][27].skip = dv28 < SKIP_VOLUME ? true : false;
      prism[i][j][28].skip = dv29 < SKIP_VOLUME ? true : false;
      prism[i][j][29].skip = dv30 < SKIP_VOLUME ? true : false;
      prism[i][j][30].skip = dv31 < SKIP_VOLUME ? true : false;
      prism[i][j][31].skip = dv32 < SKIP_VOLUME ? true : false;

      /* Calculate mass for tetrahedrons 1 to 4 */
      prism[i][j][0].dm  = dv1  * rho_th1;
      prism[i][j][1].dm  = dv2  * rho_th2;
      prism[i][j][2].dm  = dv3  * rho_th3;
      prism[i][j][3].dm  = dv4  * rho_th4;

      /* Calculate mass for tetrahedrons 5 to 8 */
      prism[i][j][4].dm  = dv5  * rho_th5;
      prism[i][j][5].dm  = dv6  * rho_th6;
      prism[i][j][6].dm  = dv7  * rho_th7;
      prism[i][j][7].dm  = dv8  * rho_th8;

      /* Calculate mass for tetrahedrons 9 to 16 */
      prism[i][j][8].dm  = dv9  * rho_th9;
      prism[i][j][9].dm  = dv10 * rho_th10;
      prism[i][j][10].dm = dv11 * rho_th11;
      prism[i][j][11].dm = dv12 * rho_th12;
      prism[i][j][12].dm = dv13 * rho_th13;
      prism[i][j][13].dm = dv14 * rho_th14;
      prism[i][j][14].dm = dv15 * rho_th15;
      prism[i][j][15].dm = dv16 * rho_th16;

      /* Calculate mass for tetrahedrons 17 to 24 */
      prism[i][j][16].dm = dv17 * rho_th17;
      prism[i][j][17].dm = dv18 * rho_th18;
      prism[i][j][18].dm = dv19 * rho_th19;
      prism[i][j][19].dm = dv20 * rho_th20;
      prism[i][j][20].dm = dv21 * rho_th21;
      prism[i][j][21].dm = dv22 * rho_th22;
      prism[i][j][22].dm = dv23 * rho_th23;
      prism[i][j][23].dm = dv24 * rho_th24;

      /* Calculate mass for tetrahedrons 25 to 32 */
      prism[i][j][24].dm = dv25 * rho_th25;
      prism[i][j][25].dm = dv26 * rho_th26;
      prism[i][j][26].dm = dv27 * rho_th27;
      prism[i][j][27].dm = dv28 * rho_th28;
      prism[i][j][28].dm = dv29 * rho_th29;
      prism[i][j][29].dm = dv30 * rho_th30;
      prism[i][j][30].dm = dv31 * rho_th31;
      prism[i][j][31].dm = dv32 * rho_th32;

      /* Assign center of mass for tetrahedrons 1 to 4 */
      prism[i][j][0].cm = cm1;
      prism[i][j][1].cm = cm2;
      prism[i][j][2].cm = cm3;
      prism[i][j][3].cm = cm4;

      /* Assign center of mass for tetrahedrons 5 to 8 */
      prism[i][j][4].cm = cm5;
      prism[i][j][5].cm = cm6;
      prism[i][j][6].cm = cm7;
      prism[i][j][7].cm = cm8;

      /* Assign center of mass for tetrahedrons 9 to 16 */
      prism[i][j][8].cm  = cm9;
      prism[i][j][9].cm  = cm10;
      prism[i][j][10].cm = cm11;
      prism[i][j][11].cm = cm12;
      prism[i][j][12].cm = cm13;
      prism[i][j][13].cm = cm14;
      prism[i][j][14].cm = cm15;
      prism[i][j][15].cm = cm16;

      /* Assign center of mass for tetrahedrons 17 to 24 */
      prism[i][j][16].cm = cm17;
      prism[i][j][17].cm = cm18;
      prism[i][j][18].cm = cm19;
      prism[i][j][19].cm = cm20;
      prism[i][j][20].cm = cm21;
      prism[i][j][21].cm = cm22;
      prism[i][j][22].cm = cm23;
      prism[i][j][23].cm = cm24;

      /* Assign center of mass for tetrahedrons 25 to 32 */
      prism[i][j][24].cm = cm25;
      prism[i][j][25].cm = cm26;
      prism[i][j][26].cm = cm27;
      prism[i][j][27].cm = cm28;
      prism[i][j][28].cm = cm29;
      prism[i][j][29].cm = cm30;
      prism[i][j][30].cm = cm31;
      prism[i][j][31].cm = cm32;
    }
}

static void sumStruct (void *in, void *inout, int *len, MPI_Datatype *type)
{
  /* Sums Triplet structures across processes using MPI reduction. */
  (void) type;

  struct Triplet *invals    = in;
  struct Triplet *inoutvals = inout;

  for (int i = 0; i < *len; i++)
  {
    /* Add x, y, and z components from input to output */
    inoutvals[i].x += invals[i].x;
    inoutvals[i].y += invals[i].y;
    inoutvals[i].z += invals[i].z;
  }
}

static void computeGravity (int ic,
                            int nl, int nf_l,
                            struct Prism prism[nl][nf_l][N_TETRAHEDRA],
                            int np, int nt,
                            struct Triplet To[np][nt], struct Triplet Go[np][nt],
                            double *mo, double GoN[np][nt])
{
  /* Computes the gravitational field and total mass using tetrahedra contributions. */
  double ml = 0; struct Triplet Gl[np][nt];

  /* Initialize arrays to zero */
  initiallizeArray (np, nt, Gl);

  int c = 0;

  clock_t starttime = clock ();

  for (int i = 0; i < nl; i++)

    for (int j = 0; j < nf_l; j++)
    {
      for (int k = 0; k < N_TETRAHEDRA; k++)
      {
        /* Skip empty tetrahedrons */
        if (prism[i][j][k].skip) continue;

        /* Get mass of tetrahedron */
        double dm = prism[i][j][k].dm;

        /* Skip tetrahedrons with zero mass */
        if (dm < EPSILON) continue;

        /* Get center of mass of tetrahedron */
        struct Triplet p, q = prism[i][j][k].cm, r;

        for (int m = 0; m < np; m++)

          for (int n = 0; n < nt; n++)
          {
            /* Set observation point for gravity calculation */
            p = To[m][n];

            /* Compute the normalized distance vector between points */
            computeDistanceVector (&p, &q, &r);

            /* Accumulate gravitational contributions in x, y, and z directions */
            Gl[m][n].x += dm * r.x;
            Gl[m][n].y += dm * r.y;
            Gl[m][n].z += dm * r.z;
          }
        
        /* Add mass of tetrahedron to total mass */
        ml += dm;
      }

      if (ic == 0) progressBar (c++, 1000, nl * nf_l, starttime);
    }

  /* Aggregate the local mass across all processes */
  MPI_Allreduce (&ml, mo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* Define an array of MPI data types for the Triplet structure */
  MPI_Datatype MPI_STRUCT, types[] = {MPI_DOUBLE,
                                      MPI_DOUBLE,
                                      MPI_DOUBLE};

  /* Create a custom MPI operation for summing Triplet structures */
  MPI_Op MPI_SUM_STRUCT;
  MPI_Op_create (sumStruct, 1, &MPI_SUM_STRUCT);

  /* Specify the number of blocks for each component of the Triplet */
  int blocks[] = {1, 1, 1};

  /* Define the byte offsets for x, y, and z components in the Triplet structure */
  MPI_Aint displ[] = {offsetof (struct Triplet, x),
                      offsetof (struct Triplet, y),
                      offsetof (struct Triplet, z)};

  /* Create a structured MPI data type for the Triplet */
  MPI_Type_create_struct (3, blocks, displ, types, &MPI_STRUCT);
  /* Commit the MPI data type for use */
  MPI_Type_commit (&MPI_STRUCT);

  struct Triplet Gt[np][nt];

  /* Reduce the gravitational contributions across all processes */
  MPI_Allreduce (Gl, Gt, np * nt, MPI_STRUCT, MPI_SUM_STRUCT, MPI_COMM_WORLD);

  /* Free the custom MPI data type */
  MPI_Type_free (&MPI_STRUCT);
  /* Free the custom MPI operation */
  MPI_Op_free (&MPI_SUM_STRUCT);

  struct Triplet Tn[np][nt];

  /* Generate normal vectors for the grid points */
  createNormalVectors (np, nt, Tn);

  for (int m = 0; m < np; m++)

    for (int n = 0; n < nt; n++)
    {
      /* Apply gravitational constant to x, y, and z components */
      Go[m][n].x = G * Gt[m][n].x;
      Go[m][n].y = G * Gt[m][n].y;
      Go[m][n].z = G * Gt[m][n].z;

      /* Compute the normal component of gravity */
      GoN[m][n] = dot (&Go[m][n], &Tn[m][n]);
    }
}

static void helpMenu (void)
{
  /* Displays the help menu with usage, arguments, description, and output files information. */
  char *help_menu = "\n GRAVITY"

                    "\n\n USAGE"
                    "\n    mpiexec -n N ./gravity DH HEIGHT"

                    "\n\n EXAMPLE"
                    "\n    mpiexec -n 4 ./gravity 1.0 255"

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    DH                    - grid spacing in degrees"
                    "\n    HEIGHT                - observation height in kilometers"

                    "\n\n DESCRIPTION"
                    "\n    Computes the Earth's gravitational field using a spherical mesh and velocity model from"
                    "\n    'input.model'. Accounts for surface topography and observation height. Reads mesh and model"
                    "\n    data from binary files (facets.bin, vertices.bin, r.bin, rho.bin, vp.bin, vs.bin) and"
                    "\n    outputs gravity components (x, y, z) and normal component to text files (G_x.dat,"
                    "\n    G_y.dat, G_z.dat, G_normal.dat). Uses MPI for parallel computation.\n\n";

  fprintf(stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  /* This function executes the gravitational field computation based on
     command-line arguments using MPI. */
  int ic;
  int nc;

  /* Initialize MPI environment */
  MPI_Init (NULL, NULL);

  /* Get the rank of the current process */
  MPI_Comm_rank (MPI_COMM_WORLD, &ic);
  /* Get the total number of processes */
  MPI_Comm_size (MPI_COMM_WORLD, &nc);

  /* Check for correct number of command-line arguments on root process */
  if (argc < 3 && ic == 0)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command line!\n");
    helpMenu ();

    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  /* Convert command-line arguments to doubles */
  double dh = atof (argv[1]);
  double height = atof (argv[2]);

  /* Calculate the number of points in longitude and latitude */
  int np = (int) (360 / dh + 1.5);
  int nt = (int) (180 / dh + 1.5);

  MPI_Barrier (MPI_COMM_WORLD);

  int ns = 0;

  /* Read number of layers and handle errors */
  if (checkMeshModelIO (readNumberOfLayers (&ns)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  int ref = 0, nf = 0, nv = 0;
  
  struct Facet **fct = NULL;
  struct Vertex **vtx = NULL;

  /* Read mesh files and handle errors */
  if (checkMeshIO (readMeshFiles (&ref, &nf, &fct, &nv, &vtx)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  /* Calculate refinement level */
  int refinement = ref + 1;

  /* Print mesh details on root process */
  if (ic == 0) fprintf (stderr, "Mesh refinement: %d\n", refinement);
  if (ic == 0) fprintf (stderr, "Number of facets: %d\n", nf);
  if (ic == 0) fprintf (stderr, "Number of vertices: %d\n", nv);

  int nf_l, shift;

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print partitioning message on root process */
  if (ic == 0) fprintf (stderr, "Partitioning facets...\n");

  partitionFacets (ic, nc, nf, &nf_l, &shift);

  struct Triplet To[np][nt];

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print grid creation message on root process */
  if (ic == 0) fprintf (stderr, "\nCreating output grid...\n");

  createOutputGrid (height, np, nt, To);

  double r[nv][ns];
  double rho[nv][ns];
  double vp[nv][ns];
  double vs[nv][ns];

  /* Print reading message on root process */
  if (ic == 0) fprintf (stderr, "Reading mesh binary files...\n");

  /* Read model files and handle errors */
  if (checkMeshModelIO (readModelFiles (nv, ns, r, rho, vp, vs)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  int nl = ns - 1;

  struct Element elm[nl][nf_l];

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print mesh creation message on root process */
  if (ic == 0) fprintf (stderr, "Creating mesh...\n");

  createMesh (nv, ns, r, rho, nf, nf_l, shift, fct[ref], vtx[ref], nl, elm);

  struct Prism prism[nl][nf_l][N_TETRAHEDRA];

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print decomposition message on root process */
  if (ic == 0) fprintf (stderr, "Decomposing prismatic elements into tetrahedra...\n");

  decomposeElements (nl, nf_l, elm, prism);

  struct Triplet Go[np][nt];

  double mo, GoN[np][nt];

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print gravity computation message on root process */
  if (ic == 0) fprintf (stderr, "Computing gravity...\n\n");

  computeGravity (ic, nl, nf_l, prism, np, nt, To, Go, &mo, GoN);

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print Earth's mass on root process */
  if (ic == 0) fprintf (stderr, "\n\nEarth's mass: %g kg\n", mo);

  MPI_Barrier (MPI_COMM_WORLD);

  /* Print writing message on root process */
  if (ic == 0) fprintf (stderr, "Writing out surfaces...\n");

  /* Write gravity data and handle errors */
  if (checkSurfaceIO (writeGravity (height, np, nt, Go, GoN)))

    MPI_Abort (MPI_COMM_WORLD, 1);

  /* Free memory allocated for facet and vertex arrays */
  destroyFacet (refinement, fct);
  destroyVertex (refinement, vtx);

  /* Print completion message on root process */
  if (ic == 0) fprintf (stderr, "Done!\n");

  MPI_Finalize ();

  return 0;
}
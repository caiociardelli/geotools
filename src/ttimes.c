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

 TTIMES

 USAGE
   ./ttimes INPUT_DIR

 EXAMPLE
   ./ttimes phases

 COMMAND-LINE ARGUMENTS
   INPUT_DIR                  - input directory containing phases_list.txt, the coordinates of all
                                ray paths for the travel times computations, and the receiver elevation

 DESCRIPTION
   Computes seismic ray travel times for phases listed in INPUT_DIR/phases_list.txt using a
   spherical mesh and velocity model. Optionally accounts for surface topography and true
   receiver elevation (read from INPUT_DIR/receiver_elevation.txt) when computing travel times.
   Reads mesh and model data from binary files (facets.bin, vertices.bin, r.bin, vp.bin, vs.bin)
   and outputs phase names and their corresponding travel times to stdout.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "structs.h"
#include "io.h"
#include "exmath.h"
#include "coordinates.h"
#include "tesselation.h"

static int resampledSize (int nd,
                          struct Triplet iray[nd],
                          double delta)
{
  /* Calculates the size of a resampled ray path based on a given spacing. */
  if (nd <= 0) return 0;
  if (delta <= 0.0) delta = TOLERANCE;

  /* Maps original indices to unique points */
  int map[nd];

  int nc = 1;
  map[0] = 0;

  for (int i = 1, j = 0; i < nd; i++)
  {
    /* Mark redundant points with same index or increment unique point count */
    if (checkRedundancy (&iray[i], &iray[i - 1]))
      
      map[i] = j;
    
      else
    {
      j++; map[i] = j; nc++;
    }
  }

  if (nc == 1) return nd;

  /* Array of unique ray points */
  struct Triplet cr[nc];
  /* Cumulative distances along ray path */
  double s[nc];

  int j = 0; cr[0] = iray[0]; s[0] = 0.0;

  for (int i = 1; i < nd; i++)
  {
    if (map[i] == j) continue;

    /* Store unique point */
    cr[++j] = iray[i];
    /* Add segment distance */
    s[j] = s[j - 1] + tripletDistance (&cr[j], &cr[j - 1]);
  }

  /* Initial size of resampled path */
  int nnd = 1;

  for (int i = 0; i + 1 < nd; i++)
  {
    /* Segment length in unique points */
    double dseg = s[map[i + 1]] - s[map[i]];

    int m = 0;

    if (dseg > delta + EPSILON)
    {
      /* Number of subdivisions needed */
      m = (int) floor ((dseg - EPSILON) / delta);
      
      if (m < 0) m = 0;
    }

    /* Add subdivisions plus endpoint to size */
    nnd += m + 1;
  }

  return nnd;
}

static void interpolateRay (int nd, struct Triplet iray[nd],
                            char ibranch[nd],
                            double delta, int nnd,
                            struct Triplet oray[nnd],
                            char obranch[nnd])
{
  /* Interpolates a ray path to achieve uniform spacing between points. */
  if (nnd <= 0) return;

  if (nd <= 0)
  {
    /* Default to P-wave for empty input */
    for (int k = 0; k < nnd; k++)
    {
      oray[k].x = 0.0;
      oray[k].y = 0.0;
      oray[k].z = 0.0;

      obranch[k] = 'P';
    }
    
    return;
  }

  /* Ensure valid resampling spacing */
  if (delta <= 0.0) delta = TOLERANCE;

  /* Maps original indices to unique points */
  int map[nd];

  /* Initialize unique point count */
  int nc = 1; map[0] = 0;

  for (int i = 1, j = 0; i < nd; i++)
  {
    /* Mark redundant points with same index or increment unique point count */
    if (checkRedundancy (&iray[i], &iray[i - 1]))
    
      map[i] = j;
    
    else
    {
      j++; map[i] = j; nc++;
    }
  }

  if (nc == 1)
  {
    /* Copy original points and branches */
    for (int i = 0; i < nd; i++)
    {
      oray[i] = iray[i];
      obranch[i] = ibranch[i];
    }

    /* Pad with last point and branch */
    for (int i = nd; i < nnd; i++)
    {
      oray[i] = iray[nd - 1];
      obranch[i] = ibranch[nd - 1];
    }

    return;
  }

  /* Array of unique ray points */
  struct Triplet cr[nc];
  /* Cumulative distances along ray path */
  double s[nc];

  cr[0] = iray[0]; s[0] = 0.0;

  for (int i = 1, j = 0; i < nd; i++)
  {
    if (map[i] == j) continue;

    /* Store unique point */
    cr[++j] = iray[i];
    /* Add segment distance */
    s[j] = s[j - 1] + tripletDistance (&cr[j], &cr[j - 1]);
  }

  /* Index for output ray points */
  int k = 0;

  /* Set first output point and branch */
  oray[k] = iray[0];
  obranch[k] = ibranch[0];
  
  k++;

  for (int i = 0; i + 1 < nd; i++)
  {
    /* Segment length in unique points */
    double si = s[map[i]];
    double sj = s[map[i + 1]];
    double dseg = sj - si;

    int m = 0;

    if (dseg > delta + EPSILON)
    
    {
      /* Number of subdivisions needed */
      m = (int) floor ((dseg - EPSILON) / delta);
    
      if (m < 0) m = 0;
    }

    for (int c = 1; c <= m; c++)
    {
      /* Target distance for interpolation */
      double su = si + c * delta;
      /* Interpolation parameter */
      double tau = 0.0;

      if (dseg > EPSILON) tau = (su - si) / dseg;

      struct Triplet p;

      struct Triplet a = cr[map[i]];
      struct Triplet b = cr[map[i + 1]];

      /* Interpolate point coordinates */
      p.x = (1.0 - tau) * a.x + tau * b.x;
      p.y = (1.0 - tau) * a.y + tau * b.y;
      p.z = (1.0 - tau) * a.z + tau * b.z;

      oray[k] = p;
      /* Assign branch type for interpolated point */
      obranch[k] = ibranch[i];
      
      k++;
    }

    /* Add endpoint of segment */
    oray[k] = iray[i + 1];
    obranch[k] = ibranch[i + 1];
    
    k++;
  }

  if (k != nnd)
  {
    /* Trim excess points if any */
    while (k > nnd) k--;
    /* Pad with last point and branch */
    while (k < nnd)
    {
      oray[k] = oray[k - 1];
      obranch[k] = obranch[k - 1];
      
      k++;
    }
  }
}

static void processPhases (char phases_path[MAX_PATH_LEN],
                           int refinement,
                           int np,
                           char phases_list[np][MAX_STRING_LEN],
                           int ref,
                           int nf, struct Facet **fct,
                           int nv, struct Vertex **vtx,
                           int ns,
                           double r[nv][ns],
                           double vp[nv][ns],
                           double vs[nv][ns],
                           double receiver_elevation)
{
  /* Computes travel times for seismic phases using a mesh and velocity model. */
  char phase_file_path[MAX_PATH_LEN];

  /* Set receiver radius */
  double receiver_radius = EARTH_RADIUS + receiver_elevation * 1E-3;
  double receiver_correction = 0.0;

  /* Set path for 1D Moho radius file */
  char moho_radius_path[MAX_PATH_LEN];

  snprintf (moho_radius_path, MAX_PATH_LEN, "input/1d_moho_radius.dat");

  /* Read 1D Moho radius and handle errors */
  int moho_index; double moho_radius;

  if (checkMohoRadiusIO (readMohoRadius (moho_radius_path, &moho_index, &moho_radius)))

    exit (EXIT_FAILURE);

  /* Process each phase in the list */
  for (int i = 0; i < np; i++)
  {
    fprintf (stderr, "phase %d: %s\n", i + 1, phases_list[i]);
    
    /* Construct path to the phase-specific ray file */
    sprintf (phase_file_path, "%s/%s.txt", phases_path, phases_list[i]);

    int nd;

    /* Read the ray path header and handle errors */
    if (checkRayIO (readRayHeader (phase_file_path, &nd)))
    
      exit (EXIT_FAILURE);

    struct Triplet iray[nd]; char ibranch[nd];

    /* Read the ray path data and handle errors */
    if (checkRayIO (readRay (phase_file_path, nd, iray, ibranch)))
    
      exit (EXIT_FAILURE);

    /* Adjust ray path for surface radius if required */
    if (INCORPORATE_SURFACE_TOPOGRAPHY)

      for (int j = 0; j < nd; j++)
      {
        struct Triplet p = iray[j];
        
        /* Calculate the current radius of the ray point */
        double point_current_radius = norm (&p);

        /* Skip points significantly deviating from Earth's surface radius */
        if (fabs (point_current_radius - EARTH_RADIUS) > THRESHOLD) continue;

        /* Skip points below the maximum Moho radius */
        if (point_current_radius < R_MOHO_MAX) continue;

        /* Normalize the point coordinates by Earth's radius */
        p.x /= EARTH_RADIUS;
        p.y /= EARTH_RADIUS;
        p.z /= EARTH_RADIUS;

        /* Determine the facet containing the current point */
        int index = findFacet (refinement, p, fct, vtx);

        /* Extract vertex indices of the containing facet */
        int i1 = fct[ref][index].i1;
        int i2 = fct[ref][index].i2;
        int i3 = fct[ref][index].i3;

        /* Retrieve radius values at the surface for each vertex */
        double r_1 = r[i1][0];
        double r_2 = r[i2][0];
        double r_3 = r[i3][0];

        /* Get the 3D coordinates of the facet vertices */
        struct Triplet p1 = vtx[ref][i1].p;
        struct Triplet p2 = vtx[ref][i2].p;
        struct Triplet p3 = vtx[ref][i3].p;

        /* Compute barycentric coordinates for the point within the facet */
        double u, v, w; baricentric (&p, &p1, &p2, &p3, &u, &v, &w);

        /* Calculate the new radius based on barycentric interpolation */
        double point_new_radius = u * r_1 + v * r_2 + w * r_3;
        /* Compute ratio of new radius to current radius */
        double ratio = point_new_radius / point_current_radius;

        /* Compute correction intervals */
        double rmin = R_MOHO_MAX;
        double rref = EARTH_RADIUS;

        /* Apply scaling factor to adjust the ray point coordinates */
        double gamma = (point_current_radius - rmin) / (rref - rmin);
        double scaling_factor = 1.0 + gamma * ratio;

        iray[j].x *= scaling_factor;
        iray[j].y *= scaling_factor;
        iray[j].z *= scaling_factor;

        /* Update receiver correction for the last point */
        if (j == nd - 1)
        
          receiver_correction = receiver_radius - point_new_radius;
      }

    /* Adjust ray path for Moho topography if required */
    if (INCORPORATE_MOHO_TOPOGRAPHY)

      for (int j = 0; j < nd; j++)
      {
        struct Triplet p = iray[j];
        
        /* Calculate the current radius of the ray point */
        double point_current_radius = norm (&p);

        /* Skip points above R_MOHO_MAX or below R_MOHO_MIN */
        if (point_current_radius > R_MOHO_MAX ||
            point_current_radius < R_MOHO_MIN) continue;

        /* Normalize the point coordinates by the current radius */
        p.x /= point_current_radius;
        p.y /= point_current_radius;
        p.z /= point_current_radius;

        /* Determine the facet containing the current point */
        int index = findFacet (refinement, p, fct, vtx);

        /* Extract vertex indices of the containing facet */
        int i1 = fct[ref][index].i1;
        int i2 = fct[ref][index].i2;
        int i3 = fct[ref][index].i3;

        /* Retrieve radius values at the Moho for each vertex */
        double r_1 = r[i1][moho_index];
        double r_2 = r[i2][moho_index];
        double r_3 = r[i3][moho_index];

        /* Get the 3D coordinates of the facet vertices */
        struct Triplet p1 = vtx[ref][i1].p;
        struct Triplet p2 = vtx[ref][i2].p;
        struct Triplet p3 = vtx[ref][i3].p;

        /* Compute barycentric coordinates for the point within the facet */
        double u, v, w; baricentric (&p, &p1, &p2, &p3, &u, &v, &w);

        /* Calculate the new moho radius based on barycentric interpolation */
        double moho_new_radius = u * r_1 + v * r_2 + w * r_3;

        /* Compute old and new Moho depths */
        double moho_depth     = EARTH_RADIUS - moho_radius;
        double moho_new_depth = EARTH_RADIUS - moho_new_radius;
        
        /* Compute ratio of Moho depth change to Moho radius */
        double ratio = (moho_depth - moho_new_depth) / moho_radius;

        /* Compute correction intervals */
        double rmin = R_MOHO_MIN;
        double rref = moho_radius;
        double rmax = R_MOHO_MAX;

        double r_point = point_current_radius;
        /* Compute scaling factor to adjust the radius */
        double gamma = (r_point > rref) ? (rmax - r_point) / (rmax - rref)
                                        : (r_point - rmin) / (rref - rmin);

        double scaling_factor = 1.0 + gamma * ratio;

        /* Apply scaling factor to adjust the ray point coordinates */
        iray[j].x *= scaling_factor;
        iray[j].y *= scaling_factor;
        iray[j].z *= scaling_factor;
      }

    int nnd = resampledSize (nd, iray, DELTA);

    struct Triplet oray[nnd]; char obranch[nnd];

    /* Interpolate ray path to uniform spacing */
    interpolateRay (nd, iray, ibranch, DELTA, nnd, oray, obranch);

    int ii[nnd];

    /* Find initial facet for the ray */
    ii[0] = findFacet (refinement, oray[0], fct, vtx);

    /* Trace ray path through facets */
    for (int i = 1; i < nnd; i++)
    {
      ii[i] = findNextFacet (ii[i - 1], &oray[i], nf, fct[ref], nv, vtx[ref]);
      
      if (ii[i] < 0)

        ii[i] = findFacet (refinement, oray[i], fct, vtx);
    }

    int jj[nnd]; jj[0] = 1;

    /* Determine shell index for each ray point */
    for (int j = 1; j < ns; j++)
    {
      int iii = ii[0];
      
      /* Extract vertex indices for the initial facet */
      int i1 = fct[ref][iii].i1;
      int i2 = fct[ref][iii].i2;
      int i3 = fct[ref][iii].i3;      

      /* Retrieve radius values for the previous and current shells */
      double r1 = r[i1][j - 1];
      double r2 = r[i2][j - 1];
      double r3 = r[i3][j - 1];
      double r4 = r[i1][j];
      double r5 = r[i2][j];
      double r6 = r[i3][j];

      /* Check if the initial ray point is inside the current element */
      if (isInsideElement (ii[0], &oray[0],
                           r1, r2, r3, r4, r5, r6,
                           nf, fct[ref], nv, vtx[ref]))
      {
        jj[0] = j; break;
      }
    }

    double inner_most_radius = r[0][ns - 2];

    /* Determine shell indices for subsequent ray points */
    for (int i = 1; i < nnd; i++)
    {
      int iii = ii[i]; jj[i] = jj[i - 1];
      int jji = jj[i - 1] - 3;
      int jjf = jj[i - 1] + 3;

      if (jji <  1) jji = 1;
      if (jjf > ns) jjf = ns;

      /* Extract vertex indices for the current facet */
      int i1 = fct[ref][iii].i1;
      int i2 = fct[ref][iii].i2;
      int i3 = fct[ref][iii].i3;

      /* Search for the appropriate shell index within a range */
      for (int j = jji; j < jjf; j++)
      {
        /* Retrieve radius values for the previous and current shells */
        double r1 = r[i1][j - 1];
        double r2 = r[i2][j - 1];
        double r3 = r[i3][j - 1];
        double r4 = r[i1][j];
        double r5 = r[i2][j];
        double r6 = r[i3][j];

        /* Calculate the radius of the current ray point */
        double radius = norm (&oray[i]);

        /* Determine the maximum radius among the previous shell vertices */
        double max_r = (r1 > r2) ? ((r1 > r3) ? r1 : r3)
                                 : ((r2 > r3) ? r2 : r3);
        /* Determine the minimum radius among the current shell vertices */
        double min_r = (r4 < r5) ? ((r4 < r6) ? r4 : r6)
                                 : ((r5 < r6) ? r5 : r6);

        if (radius < inner_most_radius)
        {
          /* Assign the innermost shell if radius is below threshold */
          jj[i] = ns - 1; break;
        }

        else if (min_r < radius &&
                 max_r > radius &&
                 isInsideElement (iii, &oray[i],
                                  r1, r2, r3, r4, r5, r6,
                                  nf, fct[ref], nv, vtx[ref]))
        {
          /* Assign the current shell index if the point is inside the element */
          jj[i] = j; break;
        }
      }
    }

    double vl[nnd];

    /* Initialize velocity at inner and outer core boundary */
    double vc1 = vp[0][ns - 1];
    double vc2 = vp[0][ns - 2];

    /* Interpolate velocity for each ray point */
    for (int i = 1; i < nnd; i++)
    {
      struct Triplet pi = oray[i - 1];
      struct Triplet pe = oray[i];

      int iii = ii[i];
      int jji = jj[i];

      /* Extract vertex indices for the current facet */
      int i1 = fct[ref][iii].i1;
      int i2 = fct[ref][iii].i2;
      int i3 = fct[ref][iii].i3;

      /* Retrieve radius values for the previous and current shells */
      double r1 = r[i1][jji - 1];
      double r2 = r[i2][jji - 1];
      double r3 = r[i3][jji - 1];
      double r4 = r[i1][jji];
      double r5 = r[i2][jji];
      double r6 = r[i3][jji];

      /* Retrieve P-wave velocity values for the previous and current shells */
      double vp1 = vp[i1][jji - 1];
      double vp2 = vp[i2][jji - 1];
      double vp3 = vp[i3][jji - 1];
      double vp4 = vp[i1][jji];
      double vp5 = vp[i2][jji];
      double vp6 = vp[i3][jji];

      /* Retrieve S-wave velocity values for the previous and current shells */
      double vs1 = vs[i1][jji - 1];
      double vs2 = vs[i2][jji - 1];
      double vs3 = vs[i3][jji - 1];
      double vs4 = vs[i1][jji];
      double vs5 = vs[i2][jji];
      double vs6 = vs[i3][jji];

      char bi = obranch[i - 1];
      char be = obranch[i];

      /* Select velocity type (P or S) based on branch at the start point */
      double vi1 = (bi == 'P') ? vp1 : vs1;
      double vi2 = (bi == 'P') ? vp2 : vs2;
      double vi3 = (bi == 'P') ? vp3 : vs3;
      double vi4 = (bi == 'P') ? vp4 : vs4;
      double vi5 = (bi == 'P') ? vp5 : vs5;
      double vi6 = (bi == 'P') ? vp6 : vs6;

      /* Select velocity type (P or S) based on branch at the end point */
      double ve1 = (be == 'P') ? vp1 : vs1; 
      double ve2 = (be == 'P') ? vp2 : vs2;
      double ve3 = (be == 'P') ? vp3 : vs3;
      double ve4 = (be == 'P') ? vp4 : vs4;
      double ve5 = (be == 'P') ? vp5 : vs5;
      double ve6 = (be == 'P') ? vp6 : vs6;

      double ri = norm (&pi);
      double re = norm (&pe);

      double vi = 0.0;
      double ve = 0.0;

      if (ri < inner_most_radius)
      {
        /* Linearly interpolate velocity for points inside the inner core */
        double m = ri / inner_most_radius;

        vi = vc1 + (vc2 - vc1) * m;
      }

      else
        /* Interpolate velocity for points outside the inner core */
        vi = interpolateElement (iii, &pi,
                                 r1, r2, r3, r4, r5, r6,
                                 vi1, vi2, vi3, vi4, vi5, vi6,
                                 nf, fct[ref], nv, vtx[ref]);

      if (re < inner_most_radius)
      {
        /* Linearly interpolate velocity for points inside the inner core */
        double m = re / inner_most_radius;

        ve = vc1 + (vc2 - vc1) * m;
      }

      else
        /* Interpolate velocity for points outside the inner core */
        ve = interpolateElement (iii, &pe,
                                 r1, r2, r3, r4, r5, r6,
                                 ve1, ve2, ve3, ve4, ve5, ve6,
                                 nf, fct[ref], nv, vtx[ref]);

      vl[i - 1] = vi; vl[i] = ve;
    }

    /* Avoid artifacts in the ray path due to mismatching velocities */
    for (int i = 0; i < nnd - 21; i++)

      if (tripletDistance (&oray[i + 10], &oray[i + 11]) < EPSILON)
      {
        for (int k = i - 10; k < i + 10; k++)

          if (fabs (vl[k] - vl[k + 1]) > 0.1)
          
            vl[k + 1] = vl[k];

        for (int k = i + 11; k < i + 20; k++)

          if (fabs (vl[k] - vl[k + 1]) > 0.1)
          
            vl[k + 1] = vl[k];
      }

    double ttime = 0.0;

    /* Calculate total travel time along the ray path */
    for (int i = 1; i < nnd; i++)
    {
      if (obranch[i - 1] != obranch[i]) continue;

      double v = 2.0 * vl[i - 1] * vl[i] / (vl[i - 1] + vl[i]);
      double d = tripletDistance (&oray[i - 1], &oray[i]);

      ttime += d / v;
    }
      
    ttime += (receiver_correction / vl[nnd - 1]);

    fprintf (stdout, "%s %.2lf\n", phases_list[i], ttime);
  }
}

static void helpMenu (void)
{
  /* Displays the help menu with usage, arguments, description, and output files information. */
  char *help_menu = "\n TTIMES"

                    "\n\n USAGE"
                    "\n    ./ttimes INPUT_DIR"

                    "\n\n EXAMPLE"
                    "\n    ./ttimes phases"

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    INPUT_DIR                  - input directory containing phases_list.txt, the coordinates of all"
                    "\n                                 ray paths for the travel times computations, and the receiver elevation"

                    "\n\n DESCRIPTION"
                    "\n    Computes seismic ray travel times for phases listed in INPUT_DIR/phases_list.txt using a"
                    "\n    spherical mesh and velocity model. Optionally accounts for surface topography and true"
                    "\n    receiver elevation (read from INPUT_DIR/receiver_elevation.txt) when computing travel times."
                    "\n    Reads mesh and model data from binary files (facets.bin, vertices.bin, r.bin, vp.bin, vs.bin)"
                    "\n    and outputs phase names and their corresponding travel times to stdout.\n\n";

  fprintf(stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  /* This function executes the travel time computation process based on command-line arguments. */
  /* Check for correct number of command-line arguments */
  if (argc < 1)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command line!\n");
    helpMenu ();

    exit (EXIT_FAILURE);
  }

  /* Initialize number of shells */
  int ns = 0;

  /* Read number of layers and handle errors */
  if (checkMeshModelIO (readNumberOfLayers (&ns)))

    exit (EXIT_FAILURE);

  /* Set path for receiver elevation file */
  char receiver_path[MAX_PATH_LEN];

  snprintf (receiver_path, MAX_PATH_LEN, "%s/receiver_elevation.txt", argv[1]);

  /* Initialize receiver elevation */
  double receiver_elevation = 0.0;

  /* Read receiver elevation if topography is incorporated */
  if (INCORPORATE_SURFACE_TOPOGRAPHY)

    if (checkReceiverIO (readReceiverElevation (&receiver_elevation, receiver_path)))
    
      exit (EXIT_FAILURE);

  int ref = 0, nf = 0, nv = 0;
  
  struct Facet **fct = NULL;
  struct Vertex **vtx = NULL;

  /* Read mesh files and handle errors */
  if (checkMeshIO (readMeshFiles (&ref, &nf, &fct, &nv, &vtx)))
  
    exit (EXIT_FAILURE);

  /* Calculate refinement level */
  int refinement = ref + 1;

  fprintf (stderr, "Mesh refinement = %d\n", refinement);
  fprintf (stderr, "Number of facets: %d\n", nf);
  fprintf (stderr, "Number of vertices: %d\n", nv);

  double r[nv][ns];
  double rho[nv][ns];
  double vp[nv][ns];
  double vs[nv][ns];

  fprintf (stderr, "Reading mesh binary files...\n");

  /* Read model files and handle errors */
  if (checkMeshModelIO (readModelFiles (nv, ns, r, rho, vp, vs)))
  
    exit (EXIT_FAILURE);

  char phases_dir[MAX_STRING_LEN];

  snprintf (phases_dir, MAX_STRING_LEN, "%s", argv[1]);

  char phases_list_path[MAX_PATH_LEN];

  snprintf (phases_list_path, MAX_PATH_LEN, "%s/phases_list.txt", phases_dir);

  int np;

  /* Read phase list header and handle errors */
  if (readPhaseHeader (phases_list_path, &np))
  {
    fprintf (stderr, "Error: could not read phases list header!\n");

    exit (EXIT_FAILURE);
  }

  char phases_list[np][MAX_STRING_LEN];

  /* Read phase list and handle errors */
  if (readPhases (phases_list_path, np, phases_list))
  {
    fprintf (stderr, "Error: could not read phases list!\n");

    exit (EXIT_FAILURE);
  }

  fprintf (stderr, "Processing phases...\n");

  /* Process all phases for travel time computation */
  processPhases (phases_dir, refinement, np, phases_list, ref,
                 nf, fct, nv, vtx, ns, r, vp, vs,
                 receiver_elevation);

  fprintf (stderr, "Done!\n");

  /* Free memory allocated for facet and vertex arrays */
  destroyFacet (refinement, fct);
  destroyVertex (refinement, vtx);

  return 0;
}
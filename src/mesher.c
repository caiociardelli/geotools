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

 MESHER

 USAGE
   ./mesher REFINEMENT MODEL

 EXAMPLE
   ./mesher 8 iasp91

 COMMAND-LINE ARGUMENTS
   REFINEMENT                 - number of mesh refinement iterations
   MODEL                      - reference 1D seismic velocity model name

 DESCRIPTION
   Creates a spherical mesh for a given 1D seismic velocity model, with optional inclusion of surface topography data
   from the ETOPO dataset at 'etopo/earth_relief_15m.topo' and 3D velocity and density lateral variations to the 1D
   reference model. The program generates a mesh based on an icosahedron with a user-specified number of refinement
   iterations, producing vertices, facets, and edges. When topography is included, it applies a 2D spherical Gaussian
   filter for anti-aliasing to adjust the topography to the mesh resolution determined by the number of refinement
   iterations. The output includes binary and text files for visualization and further analysis.

 OUTPUT FILES
   nodes.txt                   - Text file containing the longitude and latitude of mesh vertices.
   spiral.txt                  - Text file with vertex coordinates ordered to form a spiral path.
   mesh_edges_%d.txt           - Text file with coordinates of triangle edges for each refinement level.
   mesh_triangles_%d.txt       - Text file with coordinates of triangle vertices for each refinement level.
   mesh_center_triangle_%d.txt - Text file with coordinates of a single triangle's vertices.
   facets.bin                  - Binary file storing facet data (triangle indices and connectivity).
   vertices.bin                - Binary file storing vertex data (coordinates and neighbor indices).
   r.bin                       - Binary file with radius values for each vertex and shell.
   rho.bin                     - Binary file with density values for each vertex and shell.
   vp.bin                      - Binary file with P-wave velocity values for each vertex and shell.
   vs.bin                      - Binary file with S-wave velocity values for each vertex and shell.
   Mesh_VpVsRho.vtk            - VTK file for 3D visualization of the mesh with velocity and density data.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "structs.h"
#include "io.h"
#include "exmath.h"
#include "coordinates.h"
#include "tesselation.h"

/*
 * Incorporates 3D perturbations into the seismic model (radius, density, P-wave, and S-wave velocities).
 * This is a simple example demonstrating how to add 3D perturbations to the model. It applies a
 * combination of a spherical harmonic perturbation and a Gaussian pulse. For real-world applications,
 * users should implement their own function to read and interpolate an actual Earth model (e.g., from
 * a file or database) onto the mesh, replacing this example with a more accurate representation. */
static void incorparate3DModel (int nv, struct Vertex vtx[nv],
                                int ns,
                                double r[nv][ns],
                                double rho[nv][ns],
                                double vp[nv][ns],
                                double vs[nv][ns])
{
  /* Define Gaussian standard deviations for x, y, and z directions */
  double sigma_x = 0.3;
  double sigma_y = 0.3;
  double sigma_z = 0.2;

  /* Center coordinates of the Gaussian pulse */
  double x0 =  0.0;
  double y0 =  0.0;
  double z0 = -0.5;

  for (int i = 0; i < nv; i++)
  {
    double x = vtx[i].p.x;
    double y = vtx[i].p.y;
    double z = vtx[i].p.z;

    double radius, theta, phi;

    /* Convert Cartesian to spherical coordinates */
    xYZ2RThetaPhi (x, y, z, &radius, &theta, &phi);

    for (int j = 0; j < ns; j++)
    {
      /* Update radius for the current shell */
      radius = r[i][j];

      double norm_radius = radius / EARTH_RADIUS;

      /* Compute first perturbation (spherical harmonic term) */
      double dv_1 = 0.1 * (0.5 + norm_radius)
                        * sin (6 * theta) * sin (6 * phi);

      double xn = 0.0, yn = 0.0, zn = 0.0;

      /* Convert normalized radius to Cartesian coordinates */
      rThetaPhi2XYZ (norm_radius, theta, phi, &xn, &yn, &zn);

      /* Compute second perturbation (Gaussian pulse) */
      double dv_2 = -0.2 * exp (-square ((xn - x0) / sigma_x)
                                -square ((yn - y0) / sigma_y)
                                -square ((zn - z0) / sigma_z));

      double drho = (dv_1 + dv_2) * rho[i][j];
      double dvp  = (dv_1 + dv_2) * vp[i][j];
      double dvs  = (dv_1 + dv_2) * vs[i][j];

      /* Apply perturbations to the model properties */
      rho[i][j] += drho;
      vp[i][j]  += dvp;
      vs[i][j]  += dvs; 
    }
  }
}

static void filterTopography (int nlat, int nlon,
                              double lon[nlat][nlon],
                              double lat[nlat][nlon],
                              double topo[nlat][nlon],
                              int nv, struct Vertex vtx[nv],
                              int ns, double r[nv][ns],
                              double length)
{
  /* Applies a Gaussian filter to incorporate surface topography into the mesh. */
  /* Calculate Gaussian sigma based on mesh length and Earth's radius */
  double sigma = sqrt (log (2.0)) * length / (PI * EARTH_RADIUS);

  /* Process each vertex to adjust its radius based on topography */
  for (int i = 0; i < nv; i++)
  {
    double radius, theta, phi;
    /* Convert vertex coordinates to spherical */
    xYZ2RThetaPhi (vtx[i].p.x, vtx[i].p.y, vtx[i].p.z,
                   &radius, &theta, &phi);

    /* Compute center latitude in degrees */
    double lat_center = 90.0 - theta * 180.0 / PI;
    /* Compute center longitude in degrees */
    double lon_center = phi * 180.0 / PI;
    
    /* Normalize longitude to [-180, 180] */
    lon_center = fmod (lon_center + 180.0, 360.0) - 180.0;

    /* Convert sigma to degrees */
    double sigma_deg = sigma * 180.0 / PI;
    /* Define filter range as 3 times sigma in degrees */
    double delta_deg = 3.0 * sigma_deg;
    /* Define filter range in radians */
    double delta = 3.0 * sigma;

    /* Clamp minimum latitude */
    double lat_min_clamp = fmax (-90.0, lat_center - delta_deg);
    /* Clamp maximum latitude */
    double lat_max_clamp = fmin ( 90.0, lat_center + delta_deg);

    double lat0 = lat[0][0];
    /* Latitude step size */
    double dlat = lat[1][0] - lat[0][0];

    /* Starting index in float */
    double i_float_start = (lat_max_clamp - lat0) / dlat;
    /* Ending index in float */
    double i_float_end = (lat_min_clamp - lat0) / dlat;

    /* Start index */
    int i_start = fmax (0, (int) ceil (fmin (i_float_start, i_float_end)));
    /* End index */
    int i_end = fmin (nlat - 1, (int) floor (fmax (i_float_start, i_float_end)));

    if (i_start > i_end)
    {
      fprintf (stderr, "Error: invalid latitude range in function filterTopography!\n");

      exit (EXIT_FAILURE);
    }

    double lon0 = lon[0][0];
    /* Longitude step size */
    double dlon = lon[0][1] - lon[0][0];

    /* Minimum longitude of filter region */
    double lon_min = lon_center - delta_deg;
    /* Maximum longitude of filter region */
    double lon_max = lon_center + delta_deg;

    int j_start1, j_end1, j_start2 = -1, j_end2 = -1;
    
    /* Check if longitude wraps around */
    bool wraps = (lon_min < -180.0 || lon_max > 180.0);

    /* Handle longitude wrapping around the globe */
    if (lon_max - lon_min >= 360.0)
    {
      j_start1 = 0;
      j_end1 = nlon - 1;
      
      wraps = false;
    }
    
    else if (wraps)
    {
      if (lon_min < -180.0)
      {
        /* Wrapped minimum longitude */
        double lon_min_wrap = lon_min + 360.0;
        /* Wrapped maximum longitude */
        double lon_max_wrap = 180.0;
        
        j_start1 = fmax (0, (int) ceil ((lon_min_wrap - lon0) / dlon));
        j_end1 = fmin (nlon - 1, (int) floor ((lon_max_wrap - lon0) / dlon));

        /* Normal minimum longitude */
        double lon_min_normal = -180.0;
        /* Normal maximum longitude */
        double lon_max_normal = lon_max;
        
        j_start2 = fmax (0, (int) ceil ((lon_min_normal - lon0) / dlon));
        j_end2 = fmin (nlon - 1, (int) floor((lon_max_normal - lon0) / dlon));
      }
      
      else
      {
        /* Wrapped minimum longitude */
        double lon_min_wrap = -180.0;
        /* Wrapped maximum longitude */
        double lon_max_wrap = lon_max - 360.0;
        
        j_start1 = fmax (0, (int) ceil ((lon_min_wrap - lon0) / dlon));
        j_end1 = fmin (nlon - 1, (int) floor ((lon_max_wrap - lon0) / dlon));

        /* Normal minimum longitude */
        double lon_min_normal = lon_min;
        /* Normal maximum longitude */
        double lon_max_normal = 180.0;
        
        j_start2 = fmax (0, (int) ceil ((lon_min_normal - lon0) / dlon));
        j_end2 = fmin (nlon - 1, (int) floor ((lon_max_normal - lon0) / dlon));
      }
    }
    
    else
    {
      j_start1 = fmax (0, (int) ceil ((lon_min - lon0) / dlon));
      j_end1 = fmin (nlon - 1, (int) floor ((lon_max - lon0) / dlon));
    }

    /* Accumulator for weights */
    double sum_weight = 0.0;
    /* Accumulator for weighted elevations */
    double sum_topo   = 0.0;

    /* Compute weighted average elevation within the filter region */
    for (int i = i_start; i <= i_end; i++)
    {
      /* Grid point theta */
      double theta_grid = (90.0 - lat[i][0]) * PI / 180.0;

      int j_start = j_start1;
      int j_end = j_end1;
      
      for (int j = j_start; j <= j_end; j++)
      {
        /* Grid point phi */
        double phi_grid = lon[i][j] * PI / 180.0;
        /* Cosine of angular distance */
        double cos_dist = sin (theta) * sin (theta_grid) + cos (theta)
                        * cos (theta_grid) * cos (phi - phi_grid);
        
        /* Angular distance in radians */
        double dist = acos (fmax (fmin (cos_dist, 1.0), -1.0));

        if (dist > delta) continue;

        /* Gaussian weight */
        double weight = exp (-(dist * dist)
                      / (2 * square (sigma) + EPSILON));

        sum_topo   += weight * topo[i][j];
        sum_weight += weight;
      }

      if (wraps)
      {
        j_start = j_start2;
        j_end = j_end2;
        
        for (int j = j_start; j <= j_end; j++)
        {
          /* Grid point phi */
          double phi_grid = lon[i][j] * PI / 180.0;
          /* Cosine of angular distance */
          double cos_dist = sin (theta) * sin (theta_grid) + cos (theta)
                          * cos (theta_grid) * cos (phi - phi_grid);
          /* Angular distance in radians */
          double dist = acos (fmax (fmin (cos_dist, 1.0), -1.0));

          if (dist > delta) continue;

          /* Gaussian weight */
          double weight = exp (- (dist * dist)
                        / (2 * square (sigma) + EPSILON));

          sum_topo   += weight * topo[i][j];
          sum_weight += weight;
        }
      }
    }

    /* Average topographic elevation */
    double elevation = (sum_weight > 0.0) ? sum_topo / sum_weight : 0.0;
    /* Ratio of elevation to Earth radius */
    double ratio = elevation / (EARTH_RADIUS * 1E3);

    /* Compute correction intervals */
    double rmin = R_MOHO_MAX;
    double rref = EARTH_RADIUS;

    /* Apply correction */
    for (int k = 0; k < ns; k++)
    {
      if (r[i][k] < rmin) break;

      double gamma = (r[i][k] - rmin) / (rref - rmin);

      double scaling_factor = 1.0 + gamma * ratio;

      r[i][k] *= scaling_factor;

      if (r[i][k] < rmin) r[i][k] = rmin;
    }
  }
}

static void filterMoho (int nlat, int nlon,
                        double lon[nlat][nlon],
                        double lat[nlat][nlon],
                        double moho[nlat][nlon],
                        double moho_radius,
                        int nv, struct Vertex vtx[nv],
                        int ns, double r[nv][ns],
                        double length)
{
  /* Applies a Gaussian filter to incorporate surface topography into the mesh. */  
  /* Calculate Gaussian sigma based on mesh length and Earth's radius */
  double sigma = sqrt (log (2.0)) * length / (PI * moho_radius);

  /* Process each vertex to adjust its radius based on topography */
  for (int i = 0; i < nv; i++)
  {
    double radius, theta, phi;
    /* Convert vertex coordinates to spherical */
    xYZ2RThetaPhi (vtx[i].p.x, vtx[i].p.y, vtx[i].p.z,
                   &radius, &theta, &phi);

    /* Compute center latitude in degrees */
    double lat_center = 90.0 - theta * 180.0 / PI;
    /* Compute center longitude in degrees */
    double lon_center = phi * 180.0 / PI;
    
    /* Normalize longitude to [-180, 180] */
    lon_center = fmod (lon_center + 180.0, 360.0) - 180.0;

    /* Convert sigma to degrees */
    double sigma_deg = sigma * 180.0 / PI;
    /* Define filter range as 3 times sigma in degrees */
    double delta_deg = 3.0 * sigma_deg;
    /* Define filter range in radians */
    double delta = 3.0 * sigma;

    /* Clamp minimum latitude */
    double lat_min_clamp = fmax (-90.0, lat_center - delta_deg);
    /* Clamp maximum latitude */
    double lat_max_clamp = fmin ( 90.0, lat_center + delta_deg);

    double lat0 = lat[0][0];
    /* Latitude step size */
    double dlat = lat[1][0] - lat[0][0];

    /* Starting index in float */
    double i_float_start = (lat_max_clamp - lat0) / dlat;
    /* Ending index in float */
    double i_float_end = (lat_min_clamp - lat0) / dlat;

    /* Start index */
    int i_start = fmax (0, (int) ceil (fmin (i_float_start, i_float_end)));
    /* End index */
    int i_end = fmin (nlat - 1, (int) floor (fmax (i_float_start, i_float_end)));

    if (i_start > i_end)
    {
      fprintf (stderr, "Error: invalid latitude range in function filterTopography!\n");

      exit (EXIT_FAILURE);
    }

    double lon0 = lon[0][0];
    /* Longitude step size */
    double dlon = lon[0][1] - lon[0][0];

    /* Minimum longitude of filter region */
    double lon_min = lon_center - delta_deg;
    /* Maximum longitude of filter region */
    double lon_max = lon_center + delta_deg;

    int j_start1, j_end1, j_start2 = -1, j_end2 = -1;
    
    /* Check if longitude wraps around */
    bool wraps = (lon_min < -180.0 || lon_max > 180.0);

    /* Handle longitude wrapping around the globe */
    if (lon_max - lon_min >= 360.0)
    {
      j_start1 = 0;
      j_end1 = nlon - 1;
      
      wraps = false;
    }
    
    else if (wraps)
    {
      if (lon_min < -180.0)
      {
        /* Wrapped minimum longitude */
        double lon_min_wrap = lon_min + 360.0;
        /* Wrapped maximum longitude */
        double lon_max_wrap = 180.0;
        
        j_start1 = fmax (0, (int) ceil ((lon_min_wrap - lon0) / dlon));
        j_end1 = fmin (nlon - 1, (int) floor ((lon_max_wrap - lon0) / dlon));

        /* Normal minimum longitude */
        double lon_min_normal = -180.0;
        /* Normal maximum longitude */
        double lon_max_normal = lon_max;
        
        j_start2 = fmax (0, (int) ceil ((lon_min_normal - lon0) / dlon));
        j_end2 = fmin (nlon - 1, (int) floor((lon_max_normal - lon0) / dlon));
      }
      
      else
      {
        /* Wrapped minimum longitude */
        double lon_min_wrap = -180.0;
        /* Wrapped maximum longitude */
        double lon_max_wrap = lon_max - 360.0;
        
        j_start1 = fmax (0, (int) ceil ((lon_min_wrap - lon0) / dlon));
        j_end1 = fmin (nlon - 1, (int) floor ((lon_max_wrap - lon0) / dlon));

        /* Normal minimum longitude */
        double lon_min_normal = lon_min;
        /* Normal maximum longitude */
        double lon_max_normal = 180.0;
        
        j_start2 = fmax (0, (int) ceil ((lon_min_normal - lon0) / dlon));
        j_end2 = fmin (nlon - 1, (int) floor ((lon_max_normal - lon0) / dlon));
      }
    }
    
    else
    {
      j_start1 = fmax (0, (int) ceil ((lon_min - lon0) / dlon));
      j_end1 = fmin (nlon - 1, (int) floor ((lon_max - lon0) / dlon));
    }

    /* Accumulator for weights */
    double sum_weight = 0.0;
    /* Accumulator for weighted elevations */
    double sum_moho   = 0.0;

    /* Compute weighted average elevation within the filter region */
    for (int i = i_start; i <= i_end; i++)
    {
      /* Grid point theta */
      double theta_grid = (90.0 - lat[i][0]) * PI / 180.0;

      int j_start = j_start1;
      int j_end = j_end1;
      
      for (int j = j_start; j <= j_end; j++)
      {
        /* Grid point phi */
        double phi_grid = lon[i][j] * PI / 180.0;
        /* Cosine of angular distance */
        double cos_dist = sin (theta) * sin (theta_grid) + cos (theta)
                        * cos (theta_grid) * cos (phi - phi_grid);
        
        /* Angular distance in radians */
        double dist = acos (fmax (fmin (cos_dist, 1.0), -1.0));

        if (dist > delta) continue;

        /* Gaussian weight */
        double weight = exp (-(dist * dist)
                      / (2 * square (sigma) + EPSILON));

        sum_moho   += weight * moho[i][j];
        sum_weight += weight;
      }

      if (wraps)
      {
        j_start = j_start2;
        j_end = j_end2;
        
        for (int j = j_start; j <= j_end; j++)
        {
          /* Grid point phi */
          double phi_grid = lon[i][j] * PI / 180.0;
          /* Cosine of angular distance */
          double cos_dist = sin (theta) * sin (theta_grid) + cos (theta)
                          * cos (theta_grid) * cos (phi - phi_grid);
          /* Angular distance in radians */
          double dist = acos (fmax (fmin (cos_dist, 1.0), -1.0));

          if (dist > delta) continue;

          /* Gaussian weight */
          double weight = exp (- (dist * dist)
                        / (2 * square (sigma) + EPSILON));

          sum_moho   += weight * moho[i][j];
          sum_weight += weight;
        }
      }
    }

    /* Average Moho depth */
    double depth = (sum_weight > 0.0) ? sum_moho / sum_weight : 0.0;
    /* Get Moho depth for the 1D Earth model */   
    double moho_depth = EARTH_RADIUS - moho_radius;
    /* Compute ratio of Moho depth change to Moho radius */
    double ratio = (depth + moho_depth) / moho_radius;

    /* Compute correction intervals */
    double rmin = R_MOHO_MIN;
    double rref = moho_radius;
    double rmax = R_MOHO_MAX;

    /* Apply correction */
    for (int k = 0; k < ns; k++)
    {
      if (r[i][k] > rmax) continue;
      if (r[i][k] < rmin) break;

      double gamma = (r[i][k] > rref) ? (rmax - r[i][k]) / (rmax - rref)
                                      : (r[i][k] - rmin) / (rref - rmin);

      double scaling_factor = 1.0 + gamma * ratio;

      r[i][k] *= scaling_factor;

      if (r[i][k] > rmax) r[i][k] = rmax;
    }
  }
}

static void helpMenu (void)
{
  /* Displays the help menu with usage, arguments, description, and output files information. */
  char *help_menu = "\n MESHER"

                    "\n\n USAGE"
                    "\n    ./mesher REFINEMENT MODEL"

                    "\n\n EXAMPLE"
                    "\n    ./mesher 8 iasp91"

                    "\n\n COMMAND-LINE ARGUMENTS"
                    "\n    REFINEMENT                 - number of mesh refinement iterations"
                    "\n    MODEL                      - reference 1D seismic velocity model name"

                    "\n\n DESCRIPTION"
                    "\n    Creates a spherical mesh for a given 1D seismic velocity model, with optional inclusion of surface topography data"
                    "\n    from the ETOPO dataset at 'extra/earth_relief_15m.topo' and 3D velocity and density lateral variations to the 1D"
                    "\n    reference model. The program generates a mesh based on an icosahedron with a user-specified number of refinement"
                    "\n    iterations, producing vertices, facets, and edges. When topography is included, it applies a 2D spherical Gaussian"
                    "\n    filter for anti-aliasing to adjust the topography to the mesh resolution determined by the number of refinement"
                    "\n    iterations. The output includes binary and text files for visualization and further analysis."
                    "\n\n OUTPUT FILES"
                    "\n    nodes.txt                   - Text file containing the longitude and latitude of mesh vertices."
                    "\n    spiral.txt                  - Text file with vertex coordinates ordered to form a spiral path."
                    "\n    mesh_edges_%d.txt           - Text file with coordinates of triangle edges for each refinement level."
                    "\n    mesh_triangles_%d.txt       - Text file with coordinates of triangle vertices for each refinement level."
                    "\n    mesh_center_triangle_%d.txt - Text file with coordinates of a single triangle's vertices."
                    "\n    facets.bin                  - Binary file storing facet data (triangle indices and connectivity)."
                    "\n    vertices.bin                - Binary file storing vertex data (coordinates and neighbor indices)."
                    "\n    r.bin                       - Binary file with radius values for each vertex and shell."
                    "\n    rho.bin                     - Binary file with density values for each vertex and shell."
                    "\n    vp.bin                      - Binary file with P-wave velocity values for each vertex and shell."
                    "\n    vs.bin                      - Binary file with S-wave velocity values for each vertex and shell."
                    "\n    Mesh_VpVsRho.vtk            - VTK file for 3D visualization of the mesh with velocity and density data.\n\n";

  fprintf(stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  /* This function executes the meshing process based on command-line arguments. */
  /* Check for correct number of command-line arguments */
  if (argc < 3)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command line!\n");
    helpMenu ();

    exit (EXIT_FAILURE);
  }

  int refinement = atoi (argv[1]);

  /* Validate refinement value */
  if (refinement < 1)
  {
    fprintf (stderr, "Error: refinement must be greater than 0!\n");

    exit (EXIT_FAILURE);
  }

  /* Set path for topography file */
  char topo_path[MAX_PATH_LEN];

  snprintf (topo_path, MAX_PATH_LEN, "extra/earth_relief_15m.topo");

  int nlat_topo = N_LAT_ETOPO;
  int nlon_topo = N_LON_ETOPO;
  
  double lon_topo[nlat_topo][nlon_topo];
  double lat_topo[nlat_topo][nlon_topo];
  double topo[nlat_topo][nlon_topo];

  /* Read topography data and handle errors */
  if (checkTopographyIO (readTopography (nlat_topo, nlon_topo, topo_path,
                                         lon_topo, lat_topo, topo)))
    
    exit (EXIT_FAILURE);

  /* Set path for moho depth file */
  char moho_path[MAX_PATH_LEN];

  snprintf (moho_path, MAX_PATH_LEN, "extra/%s", MOHO_FILE_NAME);

  int nlat_moho = N_LAT_MOHO;
  int nlon_moho = N_LON_MOHO;
  
  double lon_moho[nlat_moho][nlon_moho];
  double lat_moho[nlat_moho][nlon_moho];
  double moho[nlat_moho][nlon_moho];

  /* Read moho depth data and handle errors */
  if (checkMohoDepthIO (readMohoDepth (nlat_moho, nlon_moho, moho_path,
                                       lon_moho, lat_moho, moho)))
    
    exit (EXIT_FAILURE);

  /* Set path for model file */
  char model_path[MAX_PATH_LEN];

  snprintf (model_path, MAX_PATH_LEN, "%s", argv[2]);

  /* Read number of shells from model header */
  int ns = 0;

  if (checkModelIO (readModelHeader (&ns, model_path)))
  
    exit (EXIT_FAILURE);

  /* Allocate memory for model arrays */
  double *rn   = malloc (ns * sizeof (double));
  double *rhon = malloc (ns * sizeof (double));
  double *vpn  = malloc (ns * sizeof (double));
  double *vsn  = malloc (ns * sizeof (double));

  if (rn == NULL || rhon == NULL || vpn == NULL || vsn == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " velocity model...\n");
    
    exit (EXIT_FAILURE);
  }

  fprintf (stderr, "Reading input model...\n");

  /* Read model data and handle errors */
  if (checkModelIO (readModel (ns, model_path,
                               rn, rhon, vpn, vsn)))
    
    exit (EXIT_FAILURE);

  fprintf (stderr, "\nCreating facets...\n");

  /* Initialize facet count */
  int nf = 0;

  /* Allocate memory for facet points */
  struct FacetPoints **fctp = createFacetPoints (refinement, &nf);

  fprintf (stderr, "Number of facets: %d\n", nf);

  /* Allocate memory for vertices and facets */
  struct Vertex **vtx = createArray (refinement);
  struct Facet  **fct = createFacets (refinement);

  /* Find nodes and facets */
  findNodes (refinement, fctp, vtx, fct);

  /* Compute number of vertices */
  int ref = refinement - 1;
  int nv = computeNumberOfVertices (ref);

  fprintf (stderr, "Number of vertices: %d\n", nv);

  /* Free memory for facet points */
  destroyFacetPoints (refinement, fctp);

  /* Find neighbors for each vertex */
  findNbVertices (refinement, fct, vtx);

  /* Write edges and facets for each refinement level if requested */
  if (DEBUG_MESH)
  {
    int index = 3; /* Index of a facet to be written */
    
    nv = 12; nf = 20;

    for (int r = 0; r < refinement; r++)
    {
      if (writeEdges (nf, fct[r], nv, vtx[r]))
      {
        fprintf (stderr, "Error: could not write edges...\n");

        exit (EXIT_FAILURE);
      }

      if (writeFacets (index, nf, fct[r], nv, vtx[r],
                      fct[r][index].nnf, fct[r][index].nbf))
      {
        fprintf (stderr, "Error: could not write facets...\n");

        exit (EXIT_FAILURE);
      }

      if (r < ref)
      {
        nf *= 4; nv = 4 * (nv - 2) + 2; index *= 4;
      }
    }
  }

  fprintf (stderr, "\nMesh statistics...\n");

  /* Compute and print mesh statistics */
  double length = meshStatistics (refinement, nf, fct[ref], nv, vtx[ref]);

  fprintf (stderr, "\nWriting out nodes...\n");

  /* Write nodes */  
  writeNodes (nv, vtx[ref]);

  fprintf (stderr, "Writing out mesh binary files...\n");  

  /* Write mesh binary files and handle errors */
  if (checkMeshModelIO (writeMeshFiles (refinement, fct, vtx)))
  
    exit (EXIT_FAILURE);

  double r[nv][ns];
  double rho[nv][ns];
  double vp[nv][ns];
  double vs[nv][ns];

  /* Initialize model arrays with 1D reference values */
  for (int i = 0; i < nv; i++)

    for (int j = 0; j < ns; j++)
    {
      r[i][j]   = rn[j];
      rho[i][j] = rhon[j];
      vp[i][j]  = vpn[j];
      vs[i][j]  = vsn[j];
    }

  /* Incorporate 3D model if specified */
  if (INCORPORATE_3D_MODEL)
  {
    fprintf (stderr, "Incorporating 3D model...\n");

    incorparate3DModel (nv, vtx[ref], ns, r, rho, vp, vs);
  }

  /* Incorporate surface topography (with anti-aliasing filter) if specified */
  if (INCORPORATE_SURFACE_TOPOGRAPHY)
  {
    fprintf (stderr, "Incorporating surface topography with anti-aliasing filter...\n");

    filterTopography (nlat_topo, nlon_topo, lon_topo, lat_topo,
                      topo, nv, vtx[ref], ns, r, length);
  }

  /* Find Moho index */
  int index = 0;

  for (int j = 1; j < ns; j++)

    if (vpn[j - 1] < 7.0 && vpn[j] > 8.0)
    {
      index = j - 1; break;
    }

  /* Get Moho radius for the 1D Reference Earth model */
  double moho_radius = r[0][index];

  /* Set path for moho radius file */
  char radius_path[MAX_PATH_LEN];

  snprintf (radius_path, MAX_PATH_LEN, "input/1d_moho_radius.dat");

  /* Write Moho radius (and associated index) for the 1D Reference Earth model
     and handle errors */
  if (checkMohoRadiusIO (writeMohoRadius (radius_path, index, moho_radius)))
  
    exit (EXIT_FAILURE);

  /* Incorporate Moho topography (with anti-aliasing filter) if specified */
  if (INCORPORATE_MOHO_TOPOGRAPHY)
  {
    fprintf (stderr, "Incorporating Moho topography with anti-aliasing filter...\n");

    filterMoho (nlat_moho, nlon_moho, lon_moho, lat_moho,
                moho, moho_radius, nv, vtx[ref], ns, r, length);
  }

  fprintf (stderr, "Writing out model binary files...\n");

  /* Write model binary files and handle errors */
  if (checkMeshModelIO (writeModelFiles (nv, ns, r, rho, vp, vs))) exit (EXIT_FAILURE);

  fprintf (stderr, "Writing out mesh in VTK format...\n");

  /* Write mesh in VTK format and handle errors */
  if (checkMeshVTK_IO (writeMeshVTK (ns, nf, fct[ref], nv, vtx[ref],
                                     r, rho, vp, vs))) exit (EXIT_FAILURE);

  fprintf (stderr, "Done!\n");

  /* Free dynamically allocated memory */
  free (rn);
  free (rhon);
  free (vpn);
  free (vsn);

  /* Free memory allocated for the facet and vertex arrays */
  destroyFacet (refinement, fct);
  destroyVertex (refinement, vtx);

  return 0;
}
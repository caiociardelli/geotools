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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "coordinates.h"

int readModelHeader (int *ns,
                     char model[MAX_STRING_LEN])
{
  /* Reads velocity model header. */
  FILE *file = fopen (model, "r");

  if (file == NULL) return 1;

  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "#N_LINES %d\n", ns) != 1) return 2;

  fclose (file);

  return 0;
}

int readModel (int ns,
               char model[MAX_STRING_LEN],
               double rn[ns], double rhon[ns],
               double vpn[ns], double vsn[ns])
{
  /* Reads radius, density, and velocities from a 1D reference model file. */
  FILE *file = fopen (model, "r");

  if (file == NULL) return 1;

  /* Skip the first three lines of the file */
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;

  double depth, rho, vp, vs;

  /* Read data for each layer */
  for (int s = 0; s < ns; s++)
  {
    if (fscanf (file, "%lf %lf %lf %lf",
                &depth, &rho, &vp, &vs) != 4) return 2;

    rn[s]   = EARTH_RADIUS - depth;
    rhon[s] = rho;
    vpn[s]  = vp;
    vsn[s]  = vs;
  }

  fclose (file);

  return 0;
}

int readTopography (int nlat, int nlon,
                    char filename[MAX_STRING_LEN],
                    double lon[nlat][nlon],
                    double lat[nlat][nlon],
                    double elev[nlat][nlon])
{
  /* Reads topography. */
  FILE *file = fopen(filename, "r");

  if (file == NULL) return 1;

  double lon_val, lat_val, elev_val;

  for (int i = 0; i < nlat; i++)
  {
    for (int j = 0; j < nlon; j++)
    {
      if (fscanf (file, "%lf %lf %lf",
                  &lon_val, &lat_val, &elev_val) != 3) return 2;

      lon[i][j]  = lon_val;
      lat[i][j]  = lat_val;
      elev[i][j] = elev_val;
    }
  }

  fclose(file);

  return 0;
}

int readReceiverElevation (double *elevation,
                           char receiver[MAX_STRING_LEN])
{
  /* Reads receiver elevation. */
  FILE *file = fopen (receiver, "r");

  if (file == NULL) return 1;

  if (fscanf (file, "#ELEVATION %lf m\n", elevation) != 1) return 2;

  fclose (file);

  return 0;
}

int readMohoRadius (char filename[MAX_STRING_LEN],
                    int *index, double *moho_radius)
{
  /* Reads moho radius for a 1D Reference Earth model. */
  FILE *file = fopen(filename, "r");

  if (file == NULL) return 1;

  /* Skip the first line of the file */
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%d %lf", index, moho_radius) != 2) return 2;

  fclose(file);

  return 0;
}

int readMohoDepth (int nlat, int nlon,
                   char filename[MAX_STRING_LEN],
                   double lon[nlat][nlon],
                   double lat[nlat][nlon],
                   double depth[nlat][nlon])
{
  /* Reads 3D moho depths. */
  FILE *file = fopen(filename, "r");

  if (file == NULL) return 1;

  double lon_val, lat_val, depth_val;

  for (int i = 0; i < nlat; i++)
  {
    for (int j = 0; j < nlon; j++)
    {
      if (fscanf (file, "%lf %lf %lf",
                  &lat_val, &lon_val, &depth_val) != 3) return 2;

      lon[i][j]   = lon_val;
      lat[i][j]   = lat_val;
      depth[i][j] = depth_val;
    }
  }

  fclose(file);

  return 0;
}

int readPhaseHeader (char *filename, int *np)
{
  /* Reads the header of a phase list file to get the number of phases. */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  char hashtag; int n;

  if (fscanf (file, "%c %d", &hashtag, &n) != 2)
  { 
    fclose (file); return 2;
  }

  fclose (file);

  *np = n;

  return 0;
}

int readPhases (char *filename,
                int np,
                char phases_list[np][MAX_STRING_LEN])
{
  /* Reads phase names from a phase list file into an array. */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  if (fscanf (file, "%*[^\n]\n") != 0) return 2;

  for (int i = 0; i < np; i++)

    if (fscanf (file, "%s", phases_list[i]) != 1)
    { 
      fclose (file); return 2;
    }

  fclose (file);

  return 0;
}

int readRayHeader (char *filename, int *nd)
{
  /* Reads the header of a ray file to get the number of ray points. */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  char hashtag; int n;

  if (fscanf (file, "%c %d", &hashtag, &n) != 2)
  { 
    fclose (file); return 2;
  }

  fclose (file);

  *nd = n;

  return 0;
}

int readRay (char *filename,
             int nd,
             struct Triplet iray[nd],
             char ibranch[nd])
{
  /* Reads ray coordinates and branch types from a ray file. */
  FILE *file = fopen (filename, "r");

  if (file == NULL) return 1;

  if (fscanf (file, "%*[^\n]\n") != 0)
  {
    fclose (file); return 2;
  }

  for (int i = 0; i < nd; i++)

    if (fscanf (file, "%lf %lf %lf %c",
                &iray[i].x, &iray[i].y, &iray[i].z,
                &ibranch[i]) != 4)
    {
      fclose (file); return 2;
    }

  fclose (file);

  return 0;
}

int readMeshFiles (int *ref,
                   int *nf, struct Facet ***fct,
                   int *nv, struct Vertex ***vtx)
{
  /* Reads mesh data from facets.bin and vertices.bin binary files. */
  char fct_name[MAX_STRING_LEN];
  char vtx_name[MAX_STRING_LEN];

  /* Set file names for facet and vertex data */
  if (sprintf (fct_name, "mesh/facets.bin") < 15) return 1;
  if (sprintf (vtx_name, "mesh/vertices.bin") < 17) return 1;

  /* Open facets binary file */
  FILE *fct_file = fopen (fct_name, "rb");

  if (fct_file == NULL) return 1;

  int refinement;

  /* Read number of refinement levels */
  if (fread (&refinement, sizeof (int), 1, fct_file) != (size_t) 1) return 2;

  *ref = refinement - 1;

  /* Allocate array for facet data across refinement levels */
  *fct = malloc (refinement * sizeof (struct Facet *));

  int read_nf = 0;

  for (int r = 0; r < refinement; r++)
  {
    int current_nf;

    /* Read number of facets for current refinement level */
    if (fread (&current_nf, sizeof (int), 1, fct_file) != (size_t) 1) return 2;

    /* Allocate memory for facets at current refinement level */
    (*fct)[r] = malloc (current_nf * sizeof (struct Facet));

    /* Read facet data for current refinement level */
    if (fread ((*fct)[r], sizeof (struct Facet), current_nf, fct_file)
               != (size_t) current_nf) return 2;

    read_nf = current_nf;
  }

  *nf = read_nf;

  /* Close facets file */
  fclose (fct_file);

  /* Open vertices binary file */
  FILE *vtx_file = fopen (vtx_name, "rb");

  if (vtx_file == NULL) return 1;

  /* Read number of refinement levels for vertices */
  if (fread (&refinement, sizeof (int), 1, vtx_file) != (size_t) 1) return 2;

  /* Check for refinement level mismatch */
  if (refinement != *ref + 1)
  {
    fprintf (stderr, "Error: refinement level mismatch between files\n");

    fclose (vtx_file);

    return 1;
  }

  /* Allocate array for vertex data across refinement levels */
  *vtx = malloc (refinement * sizeof (struct Vertex *));

  int read_nv = 0;

  for (int r = 0; r < refinement; r++)
  {
    int current_nv;

    /* Read number of vertices for current refinement level */
    if (fread (&current_nv, sizeof (int), 1, vtx_file) != (size_t) 1) return 2;

    /* Allocate memory for vertices at current refinement level */
    (*vtx)[r] = malloc (current_nv * sizeof (struct Vertex));

    /* Read vertex data for current refinement level */
    if (fread ((*vtx)[r], sizeof (struct Vertex), current_nv, vtx_file)
               != (size_t) current_nv) return 2;

    read_nv = current_nv;
  }

  *nv = read_nv;

  /* Close vertices file */
  fclose (vtx_file);

  return 0;
}

int readNumberOfLayers (int *ns)
{
  /* Reads number of layers from r.bin binary file. */
  char r_name[MAX_STRING_LEN];

  /* Set file name for data array */
  if (sprintf (r_name, "mesh/r.bin") < 10) return 1;

  /* Open radius binary file */
  FILE *r_file = fopen (r_name, "rb");
  
  if (r_file == NULL) return 2;

  int nv_r, ns_r;
  
  /* Read dimensions for radius data */
  if (fread (&nv_r, sizeof (int), 1, r_file) != (size_t) 1) return 3;
  if (fread (&ns_r, sizeof (int), 1, r_file) != (size_t) 1) return 3;

  *ns = ns_r;

  return 0;
}

int readModelFiles (int nv, int ns,
                    double r[nv][ns],
                    double rho[nv][ns],
                    double vp[nv][ns],
                    double vs[nv][ns])
{
  /* Reads model data from r.bin, rho.bin, vp.bin, and vs.bin binary files. */
  char r_name[MAX_STRING_LEN];
  char rho_name[MAX_STRING_LEN];
  char vp_name[MAX_STRING_LEN];
  char vs_name[MAX_STRING_LEN];

  /* Set file names for model data */
  if (sprintf (r_name, "mesh/r.bin") < 10) return 1;
  if (sprintf (rho_name, "mesh/rho.bin") < 12) return 1;
  if (sprintf (vp_name, "mesh/vp.bin") < 11) return 1;
  if (sprintf (vs_name, "mesh/vs.bin") < 11) return 1;

  /* Open radius binary file */
  FILE *r_file = fopen (r_name, "rb");
  
  if (r_file == NULL) return 2;
  
  int nv_r, ns_r;
  
  /* Read dimensions for radius data */
  if (fread (&nv_r, sizeof (int), 1, r_file) != (size_t) 1) return 3;
  if (fread (&ns_r, sizeof (int), 1, r_file) != (size_t) 1) return 3;

  /* Check for dimension mismatch in radius data */
  if (nv_r != nv || ns_r != ns) return 3;
  
  /* Read radius data */
  if (fread (r, sizeof (double), nv * ns, r_file) != (size_t) nv * ns) return 3;
  
  /* Close radius file */
  fclose (r_file);

  /* Open density binary file */
  FILE *rho_file = fopen (rho_name, "rb");
  
  if (rho_file == NULL) return 2;
  
  int nv_rho, ns_rho;
  
  /* Read dimensions for density data */
  if (fread (&nv_rho, sizeof (int), 1, rho_file) != (size_t) 1) return 3;
  if (fread (&ns_rho, sizeof (int), 1, rho_file) != (size_t) 1) return 3;
  
  /* Check for dimension mismatch in density data */
  if (nv_rho != nv || ns_rho != ns) return 3;
  
  /* Read density data */
  if (fread (rho, sizeof (double), nv * ns, rho_file) != (size_t) nv * ns) return 3;
  
  /* Close density file */
  fclose (rho_file);

  /* Open P-wave velocity binary file */
  FILE *vp_file = fopen (vp_name, "rb");
  
  if (vp_file == NULL) return 2;
  
  int nv_vp, ns_vp;
  
  /* Read dimensions for P-wave velocity data */
  if (fread (&nv_vp, sizeof (int), 1, vp_file) != (size_t) 1) return 3;
  if (fread (&ns_vp, sizeof (int), 1, vp_file) != (size_t) 1) return 3;
  
  /* Check for dimension mismatch in P-wave velocity data */
  if (nv_vp != nv || ns_vp != ns) return 3;
  
  /* Read P-wave velocity data */
  if (fread (vp, sizeof (double), nv * ns, vp_file) != (size_t) nv * ns) return 3;
  
  /* Close P-wave velocity file */
  fclose (vp_file);

  /* Open S-wave velocity binary file */
  FILE *vs_file = fopen (vs_name, "rb");
  
  if (vs_file == NULL) return 2;
  
  int nv_vs, ns_vs;
  
  /* Read dimensions for S-wave velocity data */
  if (fread (&nv_vs, sizeof (int), 1, vs_file) != (size_t) 1) return 3;
  if (fread (&ns_vs, sizeof (int), 1, vs_file) != (size_t) 1) return 3;
  
  /* Check for dimension mismatch in S-wave velocity data */
  if (nv_vs != nv || ns_vs != ns) return 3;
  
  /* Read S-wave velocity data */
  if (fread (vs, sizeof (double), nv * ns, vs_file) != (size_t) nv * ns) return 3;
  
  /* Close S-wave velocity file */
  fclose (vs_file);

  return 0;
}

int writeMohoRadius (char filename[MAX_STRING_LEN],
                     int index, double moho_radius)
{
  /* Writes moho radius for a 1D Reference Earth model. */
  FILE *file = fopen(filename, "w");

  if (file == NULL) return 1;

  fprintf (file, "# Index   Radius [km]\n");
  fprintf (file, "%5d %13.3lf\n", index, moho_radius);

  fclose(file);

  return 0;
}

int writeMeshFiles (int ref,
                    struct Facet **fct,
                    struct Vertex **vtx)
{
  /* Writes facet and vertex data to binary files in the mesh directory. */
  char fct_name[MAX_STRING_LEN];
  char vtx_name[MAX_STRING_LEN];

  if (sprintf (fct_name, "mesh/facets.bin") < 15) return 1;
  if (sprintf (vtx_name, "mesh/vertices.bin") < 17) return 1;

  FILE *fct_file = fopen (fct_name, "wb");

  if (fct_file == NULL) return 2;

  /* Write refinement level and facet data for each level */
  fwrite (&ref, sizeof (int), 1, fct_file);

  for (int r = 0; r < ref; r++)
  {
    int current_nf = 20 << (2 * r);

    fwrite (&current_nf, sizeof (int), 1, fct_file);
    fwrite (fct[r], sizeof (struct Facet), current_nf, fct_file);
  }

  fclose (fct_file);

  FILE *vtx_file = fopen (vtx_name, "wb");

  if (vtx_file == NULL) return 2;

  /* Write refinement level and vertex data for each level */
  fwrite (&ref, sizeof (int), 1, vtx_file);

  for (int r = 0; r < ref; r++)
  {
    int current_nv = (10 << (2 * r)) + 2;

    fwrite (&current_nv, sizeof (int), 1, vtx_file);

    fwrite (vtx[r], sizeof (struct Vertex), current_nv, vtx_file);
  }

  fclose (vtx_file);

  return 0;
}

int writeModelFiles (int nv, int ns,
                     double r[nv][ns],
                     double rho[nv][ns],
                     double vp[nv][ns],
                     double vs[nv][ns])
{
  /* Writes model data (radius, density, P-wave, and S-wave velocities) to binary files. */
  char r_name[MAX_STRING_LEN];
  char rho_name[MAX_STRING_LEN];
  char vp_name[MAX_STRING_LEN];
  char vs_name[MAX_STRING_LEN];

  if (sprintf (r_name, "mesh/r.bin") < 10) return 1;
  if (sprintf (rho_name, "mesh/rho.bin") < 12) return 1;
  if (sprintf (vp_name, "mesh/vp.bin") < 11) return 1;
  if (sprintf (vs_name, "mesh/vs.bin") < 11) return 1;

  FILE *r_file = fopen (r_name, "wb");
  
  if (r_file == NULL) return 2;
  
  /* Write radius data to file */
  fwrite (&nv, sizeof (int), 1, r_file);
  fwrite (&ns, sizeof (int), 1, r_file);
  fwrite (r, sizeof (double), nv * ns, r_file);
  
  fclose (r_file);

  FILE *rho_file = fopen (rho_name, "wb");
  
  if (rho_file == NULL) return 2;
  
  /* Write density data to file */
  fwrite (&nv, sizeof (int), 1, rho_file);
  fwrite (&ns, sizeof (int), 1, rho_file);
  fwrite (rho, sizeof (double), nv * ns, rho_file);
  
  fclose (rho_file);

  FILE *vp_file = fopen (vp_name, "wb");
  
  if (vp_file == NULL) return 2;
  
  /* Write P-wave velocity data to file */
  fwrite (&nv, sizeof (int), 1, vp_file);
  fwrite (&ns, sizeof (int), 1, vp_file);
  fwrite (vp, sizeof (double), nv * ns, vp_file);
  
  fclose (vp_file);

  FILE *vs_file = fopen (vs_name, "wb");
  
  if (vs_file == NULL) return 2;
  
  /* Write S-wave velocity data to file */
  fwrite (&nv, sizeof (int), 1, vs_file);
  fwrite (&ns, sizeof (int), 1, vs_file);
  fwrite (vs, sizeof (double), nv * ns, vs_file);
  
  fclose (vs_file);

  return 0;
}

int writeMeshVTK (int ns,
                  int nf, struct Facet fct[nf],
                  int nv, struct Vertex vtx[nv],
                  double r[nv][ns],
                  double rho[nv][ns],
                  double vp[nv][ns],
                  double vs[nv][ns])
{
  /* Writes mesh data to a VTK file for 3D visualization with velocity and density scalars. */
  char filename[MAX_STRING_LEN];

  if (sprintf (filename, "mesh/Mesh_VpVsRho.vtk") < 21) return 1;

  FILE *vtkFile = fopen (filename, "w");

  if (vtkFile == NULL) return 2;

  fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
  fprintf(vtkFile, "Mesh File\n");
  fprintf(vtkFile, "ASCII\n\n");
  fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");

  int numShells = ns;
  int numPoints = nv;

  /* Write the total number of points */
  fprintf (vtkFile, "POINTS %d float\n", numShells * numPoints);

  /* Write coordinates for all points across shells */
  for (int i = 0; i < numPoints; i++)

    for (int j = 0; j < numShells; j++)

      fprintf (vtkFile, "% .9E % .9E % .9E\n", r[i][j] * vtx[i].p.x,
                                               r[i][j] * vtx[i].p.y,
                                               r[i][j] * vtx[i].p.z);

  int numCells = nf;
  int numCellPoints = 6;

  /* Write the total number of cells and their points */
  fprintf (vtkFile, "\nCELLS %d %d\n",
           (numShells - 1) * numCells,
           (numShells - 1) * numCells * (numCellPoints + 1));

  for (int j = 0; j < numShells - 1; j++)

    for (int i = 0; i < numCells; i++)
    {
      int i1 = numShells * fct[i].i1 + j;
      int i2 = numShells * fct[i].i2 + j;
      int i3 = numShells * fct[i].i3 + j;

      fprintf (vtkFile, "%d %d %d %d %d %d %d\n", numCellPoints,
                                                  i1, i2, i3,
                                                  i1 + 1, i2 + 1, i3 + 1);
    }

  fprintf (vtkFile, "\nCELL_TYPES %d\n", (numShells - 1) * numCells);

  /* Specify cell type (13 for VTK_WEDGE) for each cell */
  for (int i = 0; i < (numShells - 1) * numCells; i++)

    fprintf (vtkFile, "13\n");

  fprintf (vtkFile, "\nPOINT_DATA %d\n", numShells * numPoints);
  fprintf (vtkFile, "SCALARS Vp float\n");
  fprintf (vtkFile, "LOOKUP_TABLE default\n");

  /* Write P-wave velocity data for each point */
  for (int i = 0; i < numPoints; i++)

    for (int j = 0; j < numShells; j++)

      fprintf (vtkFile, "%lf\n", vp[i][j]);

  fprintf (vtkFile, "SCALARS Vs float\n");
  fprintf (vtkFile, "LOOKUP_TABLE default\n");

  /* Write S-wave velocity data for each point */
  for (int i = 0; i < numPoints; i++)

    for (int j = 0; j < numShells; j++)

      fprintf (vtkFile, "%lf\n", vs[i][j]);

  fprintf (vtkFile, "SCALARS Rho float\n");
  fprintf (vtkFile, "LOOKUP_TABLE default\n");

  /* Write density data for each point */
  for (int i = 0; i < numPoints; i++)

    for (int j = 0; j < numShells; j++)

      fprintf (vtkFile, "%lf\n", rho[i][j]);

  fclose (vtkFile);

  return 0;
}

static inline void swap (int i, int j, int n,
                         int ii[n],
                         double v[n], double w[n])
{
  /* Swaps elements at indices i and j in arrays ii, v, and w. */
  int ti = ii[i]; ii[i] = ii[j]; ii[j] = ti;

  double t1 = v[i]; v[i] = v[j]; v[j] = t1;
  double t2 = w[i]; w[i] = w[j]; w[j] = t2;
}

static void insertIntoHeap (int m, int n, int ii[n],
                            double v[n], double w[n])
{
  /* Inserts a new element into the heap by comparing with its parent. */
  int p, f = m;

  while (f > 0 && v[p = f / 2] < v[f])
  {
    swap (p, f, n, ii, v, w);

    f = p;
  }
}

static void shakeHeap (int m, int n, int ii[n],
                       double v[n], double w[n])
{
  /* Adjusts the heap by moving the largest child up the tree. */
  int f = 1;

  while (f <= m)
  {
    if (f < m && v[f] < v[f + 1]) f++;
    if (v[f / 2] >= v[f]) break;

    swap (f / 2, f, n, ii, v, w);

    f *= 2;
  }
}

static void heapSort (int n, int ii[n],
                      double v[n], double w[n])
{
  /* Performs heap sort on arrays ii, v, and w using the heap data structure. */
  for (int m = 0; m < n; m++)

    insertIntoHeap (m, n, ii, v, w);

  /* Extract elements from the heap in sorted order */
  for (int m = n - 1; m > 0; m--)
  {
    swap (0, m, n, ii, v, w);
    shakeHeap (m - 1, n, ii, v, w);
  }
}

int writeNodes (int nv, struct Vertex vtx[nv])
{
  /* Writes node coordinates to nodes.txt and spiral path to spiral.txt. */
  char name1[MAX_STRING_LEN];
  char name2[MAX_STRING_LEN];

  if (sprintf (name1, "nodes.txt") < 9) return 1;
  if (sprintf (name2, "spiral.txt") < 10) return 1;

  FILE *file1 = fopen (name1, "w");
  FILE *file2 = fopen (name2, "w");

  if (file1 == NULL) return 2;
  if (file2 == NULL) return 2;

  int ii[nv];

  double tt[nv], pp[nv];

  /* Convert vertex coordinates to spherical and sort by theta */
  for (int i = 0; i < nv; i++)
  {
    ii[i] = i;

    double r, theta, phi;

    xYZ2RThetaPhi (vtx[i].p.x, vtx[i].p.y, vtx[i].p.z,
                   &r, &theta, &phi);

    if (phi < 0) phi += 2 * PI;

    tt[i] =  theta;
    pp[i] = -phi;
  }

  heapSort (nv, ii, tt, pp);

  /* Sort within theta bands by phi */
  for (int i = 0; i < nv; i++)
  {
    int j = 1;

    while (i + j < nv && fabs (tt[i] - tt[i + j]) < THRESHOLD) j++;

    heapSort (j, &ii[i], &pp[i], &tt[i]);

    i += j - 1;
  }

  /* Write sorted node coordinates to nodes.txt */
  for (int i = 0; i < nv; i++)
  {
    double r, theta, phi;

    xYZ2RThetaPhi (vtx[ii[i]].p.x, vtx[ii[i]].p.y, vtx[ii[i]].p.z,
                   &r, &theta, &phi);

    double lat, lon; capCoordinates (theta, phi, &lat, &lon);

    fprintf (file1, "%lf %lf\n", lon, lat);
  }

  /* Write spiral path coordinates to spiral.txt */
  for (int i = 1; i < nv; i++)
  {
    double r, theta, phi, lat, lon;

    xYZ2RThetaPhi (vtx[ii[i - 1]].p.x, vtx[ii[i - 1]].p.y, vtx[ii[i - 1]].p.z,
                   &r, &theta, &phi);

    capCoordinates (theta, phi, &lat, &lon);

    fprintf (file2, "%lf %lf\n", lon, lat);

    xYZ2RThetaPhi (vtx[ii[i]].p.x, vtx[ii[i]].p.y, vtx[ii[i]].p.z,
                   &r, &theta, &phi);

    capCoordinates (theta, phi, &lat, &lon);

    fprintf (file2, "%lf %lf\n>\n", lon, lat);
  }

  fclose (file1);

  return 0;
}

int writeEdges (int nf, struct Facet fct[nf],
                int nv, struct Vertex vtx[nv])
{
  /* Writes edge coordinates for all facets to a file. */
  char name[MAX_STRING_LEN];

  sprintf (name, "mesh_edges_%d.txt", nf);

  FILE *file = fopen (name, "w");

  if (file == NULL) return 1;

  /* Write coordinates for each edge of every facet */
  for (int j = 0; j < nf; j++)
  {
    double r, theta, phi;

    xYZ2RThetaPhi (vtx[fct[j].i1].p.x,
                   vtx[fct[j].i1].p.y,
                   vtx[fct[j].i1].p.z,
                   &r, &theta, &phi);

    double lat1, lon1; capCoordinates (theta, phi, &lat1, &lon1);

    xYZ2RThetaPhi (vtx[fct[j].i2].p.x,
                   vtx[fct[j].i2].p.y,
                   vtx[fct[j].i2].p.z,
                   &r, &theta, &phi);

    double lat2, lon2; capCoordinates (theta, phi, &lat2, &lon2);

    xYZ2RThetaPhi (vtx[fct[j].i3].p.x,
                   vtx[fct[j].i3].p.y,
                   vtx[fct[j].i3].p.z,
                   &r, &theta, &phi);

    double lat3, lon3; capCoordinates (theta, phi, &lat3, &lon3);

    fprintf (file, "%lf %lf\n",    lon1, lat1);
    fprintf (file, "%lf %lf\n",    lon2, lat2);
    fprintf (file, "%lf %lf\n",    lon3, lat3);
    fprintf (file, "%lf %lf\n>\n", lon1, lat1);
  }

  fclose (file);

  return 0;
}

int writeFacets (int index,
                 int nf, struct Facet fct[nf],
                 int nv, struct Vertex vtx[nv],
                 int nn, int nbf[nn])
{
  /* Writes facet coordinates to a file and a single facet to a separate file. */
  char name1[MAX_STRING_LEN];
  char name2[MAX_STRING_LEN];

  sprintf (name1, "mesh_triangles_%d.txt", nf);
  sprintf (name2, "mesh_center_triangle_%d.txt", nf);

  FILE *file1 = fopen (name1, "w");
  FILE *file2 = fopen (name2, "w");

  if (file1 == NULL) return 1;
  if (file2 == NULL) return 1;

  /* Write coordinates for all specified facets to the first file */
  for (int i = 0; i < nn; i++)
  {
    int ii = nbf[i];

    double r, theta, phi;

    xYZ2RThetaPhi (vtx[fct[ii].i1].p.x,
                   vtx[fct[ii].i1].p.y,
                   vtx[fct[ii].i1].p.z,
                   &r, &theta, &phi);

    double lat1, lon1; capCoordinates (theta, phi, &lat1, &lon1);

    xYZ2RThetaPhi (vtx[fct[ii].i2].p.x,
                   vtx[fct[ii].i2].p.y,
                   vtx[fct[ii].i2].p.z,
                   &r, &theta, &phi);

    double lat2, lon2; capCoordinates (theta, phi, &lat2, &lon2);

    xYZ2RThetaPhi (vtx[fct[ii].i3].p.x,
                   vtx[fct[ii].i3].p.y,
                   vtx[fct[ii].i3].p.z,
                   &r, &theta, &phi);

    double lat3, lon3; capCoordinates (theta, phi, &lat3, &lon3);

    fprintf (file1, "%lf %lf\n",    lon1, lat1);
    fprintf (file1, "%lf %lf\n",    lon2, lat2);
    fprintf (file1, "%lf %lf\n",    lon3, lat3);
    fprintf (file1, "%lf %lf\n>\n", lon1, lat1);
  }

  fclose (file1);

  int ii = index;

  /* Write coordinates for a single specified facet to the second file */
  double r, theta, phi;

  xYZ2RThetaPhi (vtx[fct[ii].i1].p.x,
                 vtx[fct[ii].i1].p.y,
                 vtx[fct[ii].i1].p.z,
                 &r, &theta, &phi);

  double lat1, lon1; capCoordinates (theta, phi, &lat1, &lon1);

  xYZ2RThetaPhi (vtx[fct[ii].i2].p.x,
                 vtx[fct[ii].i2].p.y,
                 vtx[fct[ii].i2].p.z,
                 &r, &theta, &phi);

  double lat2, lon2; capCoordinates (theta, phi, &lat2, &lon2);

  xYZ2RThetaPhi (vtx[fct[ii].i3].p.x,
                 vtx[fct[ii].i3].p.y,
                 vtx[fct[ii].i3].p.z,
                 &r, &theta, &phi);

  double lat3, lon3; capCoordinates (theta, phi, &lat3, &lon3);

  fprintf (file2, "%lf %lf\n",    lon1, lat1);
  fprintf (file2, "%lf %lf\n",    lon2, lat2);
  fprintf (file2, "%lf %lf\n",    lon3, lat3);
  fprintf (file2, "%lf %lf\n>\n", lon1, lat1);

  fclose (file2);

  return 0;
}

int writeGravity (double height,
                  int np, int nt,
                  struct Triplet Go[np][nt],
                  double GoN[np][nt])
{
  /* Writes gravitational field components to output files. */
  char name1[MAX_STRING_LEN];
  char name2[MAX_STRING_LEN];
  char name3[MAX_STRING_LEN];
  char name4[MAX_STRING_LEN];

  double Theta[nt], Phi[np];

  double p = -180;
  double t =   90;

  double dp =  360 / (np - 1);
  double dt = -180 / (nt - 1);

  for (int m = 0; m < np; m++)
  {
    /* Populate longitude array with incremental steps */
    Phi[m] = p; p += dp;
  }

  for (int n = 0; n < nt; n++)
  {
    /* Populate latitude array with incremental steps */
    Theta[n] = t; t += dt;
  }

  FILE *file1, *file2, *file3, *file4;

  /* Create output file names */
  if (sprintf (name1, "G_x.dat") < 7) return 1;
  if (sprintf (name2, "G_y.dat") < 7) return 1;
  if (sprintf (name3, "G_z.dat") < 7) return 1;
  if (sprintf (name4, "G_normal.dat") < 12) return 1;

  /* Open output files */
  file1 = fopen (name1, "w");
  file2 = fopen (name2, "w");
  file3 = fopen (name3, "w");
  file4 = fopen (name4, "w");

  /* Check if files were successfully opened */
  if (file1 == NULL) return 2;
  if (file2 == NULL) return 2;
  if (file3 == NULL) return 2;
  if (file4 == NULL) return 2;

  /* Write header information to output files */
  fprintf (file1, "#height (km) nlat nlon: %lg %u %u\n", height, nt, np);
  fprintf (file2, "#height (km) nlat nlon: %lg %u %u\n", height, nt, np);
  fprintf (file3, "#height (km) nlat nlon: %lg %u %u\n", height, nt, np);
  fprintf (file4, "#height (km) nlat nlon: %lg %u %u\n", height, nt, np);
  fprintf (file1, "#latitude (degrees)   longitude (degrees)     g(m/s^2)\n");
  fprintf (file2, "#latitude (degrees)   longitude (degrees)     g(m/s^2)\n");
  fprintf (file3, "#latitude (degrees)   longitude (degrees)     g(m/s^2)\n");
  fprintf (file4, "#latitude (degrees)   longitude (degrees)     g(m/s^2)\n");

  for (int n = 0; n < nt; n++)

    for (int m = 0; m < np; m++)
    {
      /* Write x, y, z, and normal components to respective files */
      fprintf (file1, "%15.7lf %21.7lf %18lE\n", Theta[n], Phi[m],
                                                 Go[m][n].x);
      fprintf (file2, "%15.7lf %21.7lf %18lE\n", Theta[n], Phi[m],
                                                 Go[m][n].y);
      fprintf (file3, "%15.7lf %21.7lf %18lE\n", Theta[n], Phi[m],
                                                 Go[m][n].z);
      fprintf (file4, "%15.7lf %21.7lf %18lE\n", Theta[n], Phi[m],
                                                 GoN[m][n]);
    }

  /* Close output files */
  fclose (file1);
  fclose (file2);
  fclose (file3);
  fclose (file4);

  return 0;
}

int checkModelIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open 1D reference model file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read 1D reference model file!\n");
    break;
  }

  return rvalue;
}

int checkTopographyIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open topography file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read topography file!\n");
    break;
  }

  return rvalue;
}

int checkReceiverIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open receiver_elevation file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read receiver elevation file!\n");
    break;
  }

  return rvalue;
}

int checkMohoRadiusIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open 1D Moho radius file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read 1D Moho radius file!\n");
    break;
  }

  return rvalue;
}

int checkMohoDepthIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open Moho depth file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read Moho depth file!\n");
    break;
  }

  return rvalue;
}

int checkMeshIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not open file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file!\n");
    break;
  }

  return rvalue;
}

int checkMeshModelIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not open file!\n");
    break;

    case 3:
      fprintf (stderr, "Error: could not read file. Dimension mismatch!\n");
    break;
  }

  return rvalue;
}

int checkMeshVTK_IO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: unable to write file name to string buffer!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not open file!\n");
    break;
  }

  return rvalue;
}

int checkRayIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open ray coordinates file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read ray coordinates file!\n");
    break;
  }

  return rvalue;
}

int checkSurfaceIO (int rvalue)
{
  /* Prints error messages for the I/O operations and returns the error code. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Unable to write file name to string buffer...\n");
    break;

    case 2:
      fprintf (stderr, "Error opening file...\n");
    break;
  }

  return rvalue;
}

double meshStatistics (int refinement,
                       int nf, struct Facet fct[nf],
                       int nv, struct Vertex vtx[nv])
{
  /* Computes and prints statistical measures (min, max, avg) of triangle side lengths in the mesh. */
  double min_arc =  INFINITY;
  double max_arc = -INFINITY;
  double avg_arc =  0.0;

  /* Calculate distances for each triangle and update statistics */
  for (int i = 0; i < nf; i++)
  {
    struct Triplet p1 = vtx[fct[i].i1].p;
    struct Triplet p2 = vtx[fct[i].i2].p;
    struct Triplet p3 = vtx[fct[i].i3].p;

    double d1 = EARTH_RADIUS * angle (&p1, &p2);
    double d2 = EARTH_RADIUS * angle (&p1, &p3);
    double d3 = EARTH_RADIUS * angle (&p2, &p3);

    if (d1 > max_arc) max_arc = d1;
    if (d2 > max_arc) max_arc = d2;
    if (d3 > max_arc) max_arc = d3;

    if (d1 < min_arc) min_arc = d1;
    if (d2 < min_arc) min_arc = d2;
    if (d3 < min_arc) min_arc = d3;

    avg_arc += d1 + d2 + d3;
  }

  /* Compute average and print all statistics */
  avg_arc /= (3 * nf);

  fprintf (stderr, "Refinement = %d\n", refinement);
  fprintf (stderr, "Number of vertices for whole sphere would be: %d\n", nv);
  fprintf (stderr, "Minimum triangle-side length, %lf, km\n", min_arc);
  fprintf (stderr, "Average triangle-side length, %lf, km\n", avg_arc);
  fprintf (stderr, "Maximum triangle-side length, %lf, km\n", max_arc);

  return avg_arc;
}
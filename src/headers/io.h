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

#ifndef IO_H
#define IO_H

#include "constants.h"
#include "structs.h"

int readModelHeader (int *ns,
                     char model[MAX_STRING_LEN]);
int readModel (int ns,
               char model[MAX_STRING_LEN],
               double rn[ns], double rhon[ns],
               double vpn[ns], double vsn[ns]);
int readTopography (int nlat, int nlon,
                    char filename[MAX_STRING_LEN],
                    double lon[nlat][nlon],
                    double lat[nlat][nlon],
                    double elev[nlat][nlon]);
int readReceiverElevation (double *elevation,
                           char receiver[MAX_STRING_LEN]);
int readPhaseHeader (char *filename, int *np);
int readPhases (char *filename,
                int np,
                char phases_list[np][MAX_STRING_LEN]);
int readRayHeader (char *filename, int *nd);
int readRay (char *filename,
             int nd,
             struct Triplet iray[nd],
             char ibranch[nd]);
int readMeshFiles (int *ref,
                   int *nf, struct Facet ***fct,
                   int *nv, struct Vertex ***vtx);
int readNumberOfLayers (int *ns);
int readModelFiles (int nv, int ns,
                    double r[nv][ns],
                    double rho[nv][ns],
                    double vp[nv][ns],
                    double vs[nv][ns]);

int writeMeshFiles (int ref,
                    struct Facet **fct,
                    struct Vertex **vtx);
int writeModelFiles (int nv, int ns,
                     double r[nv][ns],
                     double rho[nv][ns],
                     double vp[nv][ns],
                     double vs[nv][ns]);
int writeMeshVTK (int ns,
                  int nf, struct Facet fct[nf],
                  int nv, struct Vertex vtx[nv],
                  double r[nv][ns],
                  double rho[nv][ns],
                  double vp[nv][ns],
                  double vs[nv][ns]);
int writeNodes (int nv, struct Vertex vtx[nv]);
int writeEdges (int nf, struct Facet fct[nf],
                int nv, struct Vertex vtx[nv]);
int writeFacets (int index,
                 int nf, struct Facet fct[nf],
                 int nv, struct Vertex vtx[nv],
                 int nn, int nbf[nn]);
int writeGravity (double height,
                  int np, int nt,
                  struct Triplet Go[np][nt],
                  double GoN[np][nt]);

int checkModelIO (int rvalue);
int checkTopographyIO (int rvalue);
int checkReceiverIO (int rvalue);
int checkMeshIO (int rvalue);
int checkMeshModelIO (int rvalue);
int checkMeshVTK_IO (int rvalue);
int checkRayIO (int rvalue);
int checkSurfaceIO (int rvalue);

double meshStatistics (int refinement,
                       int nf, struct Facet fct[nf],
                       int nv, struct Vertex vtx[nv]);
#endif
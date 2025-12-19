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

#ifndef TESSELATION_H
#define TESSELATION_H

#include "constants.h"
#include "structs.h"

int computeNumberOfVertices (int ref);

bool notInList (int index, int n,
                int nbf[N_MAX_NGB_VERTICES]);
bool isInsideFacet (int i,
                    struct Triplet *p,
                    int nf, struct Facet fct[nf],
                    int nv, struct Vertex vtx[nv]);
bool isInsideElement (int i, struct Triplet *p,
                      double r1, double r2, double r3,
                      double r4, double r5, double r6,
                      int nf, struct Facet fct[nf],
                      int nv, struct Vertex vtx[nv]);

struct FacetPoints **createFacetPoints (int refinement, int *nf);
struct Vertex **createArray (int refinement);
struct Facet **createFacets (int refinement);

void destroyFacetPoints (int refinement, struct FacetPoints **fctp);
void destroyFacet (int refinement, struct Facet **fct);
void destroyVertex (int refinement, struct Vertex **vtx);

void findNbFacets (int i, int nn, int nbf[nn], int nf,
                   struct FacetPoints fctp[nf]);
void findNodes (int refinement,
                struct FacetPoints **fctp,
                struct Vertex **vtx,
                struct Facet **fct);
void findNbVertices (int refinement,
                     struct Facet **fct,
                     struct Vertex **vtx);

int findNextFacet (int i,
                   struct Triplet *p,
                   int nf, struct Facet fct[nf],
                   int nv, struct Vertex vtx[nv]);
int _findFacet (int refinement,
                int r, int i, int nf, int nv,
                struct Triplet p,
                struct Facet **fct,
                struct Vertex **vtx);
int findFacet (int refinement, 
               struct Triplet p,
               struct Facet **fct,
               struct Vertex **vtx);

double interpolateElement (int i, struct Triplet *p,
                           double r1, double r2, double r3,
                           double r4, double r5, double r6,
                           double m1, double m2, double m3,
                           double m4, double m5, double m6,
                           int nf, struct Facet fct[nf],
                           int nv, struct Vertex vtx[nv]);
#endif
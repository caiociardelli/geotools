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

#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "exmath.h"
#include "coordinates.h"
#include "tesselation.h"

int computeNumberOfVertices (int ref)
{
  /* Computes the total number of vertices after a given number of refinements. */
  int nv = 12;

  for (int r = 0; r < ref; r++)

    nv = 4 * (nv - 2) + 2;

  return nv;
}

bool notInList (int index, int n,
                int nbf[N_MAX_NGB_VERTICES])
{
  /* Checks if an index is not already in the list of neighboring vertices. */
  for (int i = 0; i < n; i++)

    if (index == nbf[i]) return false;

  return true;
}

bool isInsideFacet (int i,
                    struct Triplet *p,
                    int nf, struct Facet fct[nf],
                    int nv, struct Vertex vtx[nv])
{
  /* Checks if a point lies inside a triangular facet. */
  struct Triplet p1 = vtx[fct[i].i1].p;
  struct Triplet p2 = vtx[fct[i].i2].p;
  struct Triplet p3 = vtx[fct[i].i3].p;

  struct Triplet q = project (p, &p1, &p2, &p3);

  double u, v, w; baricentric (&q, &p1, &p2, &p3, &u, &v, &w);

  if (fabs (u) +
      fabs (v) +
      fabs (w) < LEEWAY)

    return true;

  return false;
}

bool isInsideElement (int i, struct Triplet *p,
                      double r1, double r2, double r3,
                      double r4, double r5, double r6,
                      int nf, struct Facet fct[nf],
                      int nv, struct Vertex vtx[nv])
{
  /* Checks if a point lies inside a prismatic element. */
  double w1, w2, w3, w4, w5, w6;

  weights (p,
           r1, r2, r3, r4, r5, r6,
           &fct[i], nv, vtx,
           &w1, &w2, &w3, &w4, &w5, &w6);

  double w = fabs (w1) + fabs (w2) + fabs (w3)
           + fabs (w4) + fabs (w5) + fabs (w6);

  if (w < LEEWAY) return true;

  return false;
}

struct FacetPoints **createFacetPoints (int refinement, int *nf)
{
  /* Creates and initializes an array of facet points for multiple refinement levels. */
  int Nv = 12, Nf = 20;

  struct FacetPoints **fctp = malloc (refinement * sizeof (struct FacetPoints*));

  /* Allocate memory for each refinement level */
  for (int i = 0, n = Nf; i < refinement; i++, n *= 4)

    fctp[i] = malloc (n * sizeof (struct FacetPoints));

  /* Compute the golden ratio for icosahedron vertex placement */
  double PHI = goldenRatio ();

  /* Array to store initial icosahedron vertex coordinates */
  struct Triplet Icv[Nv];

  Icv[0].x  =  PHI; Icv[0].y  =  1.0; Icv[0].z  =  0.0; /* Vertex 0 coordinates */
  Icv[1].x  =  PHI; Icv[1].y  = -1.0; Icv[1].z  =  0.0; /* Vertex 1 coordinates */
  Icv[2].x  = -PHI; Icv[2].y  = -1.0; Icv[2].z  =  0.0; /* Vertex 2 coordinates */
  Icv[3].x  = -PHI; Icv[3].y  =  1.0; Icv[3].z  =  0.0; /* Vertex 3 coordinates */
  Icv[4].x  =  1.0; Icv[4].y  =  0.0; Icv[4].z  =  PHI; /* Vertex 4 coordinates */
  Icv[5].x  = -1.0; Icv[5].y  =  0.0; Icv[5].z  =  PHI; /* Vertex 5 coordinates */
  Icv[6].x  = -1.0; Icv[6].y  =  0.0; Icv[6].z  = -PHI; /* Vertex 6 coordinates */
  Icv[7].x  =  1.0; Icv[7].y  =  0.0; Icv[7].z  = -PHI; /* Vertex 7 coordinates */
  Icv[8].x  =  0.0; Icv[8].y  =  PHI; Icv[8].z  =  1.0; /* Vertex 8 coordinates */
  Icv[9].x  =  0.0; Icv[9].y  =  PHI; Icv[9].z  = -1.0; /* Vertex 9 coordinates */
  Icv[10].x =  0.0; Icv[10].y = -PHI; Icv[10].z = -1.0; /* Vertex 10 coordinates */
  Icv[11].x =  0.0; Icv[11].y = -PHI; Icv[11].z =  1.0; /* Vertex 11 coordinates */

  /* Rotation angle around x-axis */
  double alpha = 0.0;
  /* Rotation angle around y-axis based on vertex angle */
  double beta  = 0.5 * angle (&Icv[4], &Icv[5]);
  /* Rotation angle around z-axis */
  double gamma = 0.0;
  
  /* Apply rotation to initial vertices */
  rotate (Nv, Icv, alpha, beta, gamma);
  
  /* Array to store initial facet indices */
  struct Indices jj[Nf];

  jj[0].i1  = 4;  jj[0].i2  = 5;  jj[0].i3  = 8;        /* Facet 0 indices */
  jj[1].i1  = 5;  jj[1].i2  = 4;  jj[1].i3  = 11;       /* Facet 1 indices */
  jj[2].i1  = 2;  jj[2].i2  = 5;  jj[2].i3  = 11;       /* Facet 2 indices */
  jj[3].i1  = 8;  jj[3].i2  = 5;  jj[3].i3  = 3;        /* Facet 3 indices */
  jj[4].i1  = 0;  jj[4].i2  = 4;  jj[4].i3  = 8;        /* Facet 4 indices */
  jj[5].i1  = 1;  jj[5].i2  = 4;  jj[5].i3  = 0;        /* Facet 5 indices */
  jj[6].i1  = 11; jj[6].i2  = 4;  jj[6].i3  = 1;        /* Facet 6 indices */
  jj[7].i1  = 10; jj[7].i2  = 11; jj[7].i3  = 1;        /* Facet 7 indices */
  jj[8].i1  = 2;  jj[8].i2  = 11; jj[8].i3  = 10;       /* Facet 8 indices */
  jj[9].i1  = 3;  jj[9].i2  = 5;  jj[9].i3  = 2;        /* Facet 9 indices */
  jj[10].i1 = 6;  jj[10].i2 = 2;  jj[10].i3 = 10;       /* Facet 10 indices */
  jj[11].i1 = 6;  jj[11].i2 = 3;  jj[11].i3 = 2;        /* Facet 11 indices */
  jj[12].i1 = 9;  jj[12].i2 = 3;  jj[12].i3 = 6;        /* Facet 12 indices */
  jj[13].i1 = 9;  jj[13].i2 = 8;  jj[13].i3 = 3;        /* Facet 13 indices */
  jj[14].i1 = 0;  jj[14].i2 = 8;  jj[14].i3 = 9;        /* Facet 14 indices */
  jj[15].i1 = 7;  jj[15].i2 = 1;  jj[15].i3 = 0;        /* Facet 15 indices */
  jj[16].i1 = 10; jj[16].i2 = 1;  jj[16].i3 = 7;        /* Facet 16 indices */
  jj[17].i1 = 7;  jj[17].i2 = 0;  jj[17].i3 = 9;        /* Facet 17 indices */
  jj[18].i1 = 6;  jj[18].i2 = 10; jj[18].i3 = 7;        /* Facet 18 indices */
  jj[19].i1 = 6;  jj[19].i2 = 7;  jj[19].i3 = 9;        /* Facet 19 indices */

  /* Assign initial facet points and compute normals */
  for (int j = 0; j < Nf; j++)
  {
    fctp[0][j].p1 = Icv[jj[j].i1];
    fctp[0][j].p2 = Icv[jj[j].i2];
    fctp[0][j].p3 = Icv[jj[j].i3];

    normalize (&fctp[0][j].p1);
    normalize (&fctp[0][j].p2);
    normalize (&fctp[0][j].p3);

    fctp[0][j].n_sph = normalVector (&fctp[0][j].p1,
                                     &fctp[0][j].p2,
                                     &fctp[0][j].p3);

    fctp[0][j].min_dot = dot (&fctp[0][j].n_sph, &fctp[0][j].p1);
  }

  int nbfi[] = {0, 1, 2, 3, 4,
                5, 6, 7, 8, 9,
                10, 11, 12, 13,
                14, 15, 16, 17,
                18, 19};

  /* Find initial neighboring facets */
  for (int n = 0; n < Nf; n++)

    findNbFacets (n, Nf, nbfi, Nf, fctp[0]);

  int nn = 4 * N_MAX_NGB_FACETS, nbf[nn];

  /* Refine facets for each level */
  for (int i = 1; i < refinement; i++)
  {
    for (int j = 0; j < Nf; j++)
    {
      int k = 4 * j;

      fctp[i][k + 0] = fctp[i - 1][j];
      fctp[i][k + 1] = fctp[i - 1][j];
      fctp[i][k + 2] = fctp[i - 1][j];
      fctp[i][k + 3] = fctp[i - 1][j];

      struct Triplet p1 = midPoint (&fctp[i - 1][j].p1, &fctp[i - 1][j].p2);
      struct Triplet p2 = midPoint (&fctp[i - 1][j].p2, &fctp[i - 1][j].p3);
      struct Triplet p3 = midPoint (&fctp[i - 1][j].p3, &fctp[i - 1][j].p1);

      normalize (&p1);
      normalize (&p2);
      normalize (&p3);

      fctp[i][k].p2     = p1;
      fctp[i][k].p3     = p3;
      fctp[i][k + 1].p1 = p1;
      fctp[i][k + 1].p3 = p2;
      fctp[i][k + 2].p1 = p3;
      fctp[i][k + 2].p2 = p2;
      fctp[i][k + 3].p1 = p1;
      fctp[i][k + 3].p2 = p2;
      fctp[i][k + 3].p3 = p3;

      fctp[i][k].n_sph = normalVector (&fctp[i][k].p1,
                                       &fctp[i][k].p2,
                                       &fctp[i][k].p3);

      fctp[i][k].min_dot = dot (&fctp[i][k].n_sph, &fctp[i][k].p1);

      fctp[i][k + 1].n_sph = normalVector (&fctp[i][k + 1].p1,
                                           &fctp[i][k + 1].p2,
                                           &fctp[i][k + 1].p3);

      fctp[i][k + 1].min_dot = dot (&fctp[i][k + 1].n_sph, &fctp[i][k + 1].p1);

      fctp[i][k + 2].n_sph = normalVector (&fctp[i][k + 2].p1,
                                           &fctp[i][k + 2].p2,
                                           &fctp[i][k + 2].p3);

      fctp[i][k + 2].min_dot = dot (&fctp[i][k + 2].n_sph, &fctp[i][k + 2].p1);

      fctp[i][k + 3].n_sph = normalVector (&fctp[i][k + 3].p1,
                                           &fctp[i][k + 3].p2,
                                           &fctp[i][k + 3].p3);

      fctp[i][k + 3].min_dot = dot (&fctp[i][k + 3].n_sph, &fctp[i][k + 3].p1);
    }

    /* Update neighboring facets for refined triangles */
    for (int j = 0; j < Nf; j++)
    {
      int k = 4 * j, nn = 4 * fctp[i - 1][j].nnf;

      for (int l = 0, m = 0; l < nn; l += 4, m++)
      {
        int n = 4 * fctp[i - 1][j].nbf[m];

        nbf[l]     = n;
        nbf[l + 1] = n + 1;
        nbf[l + 2] = n + 2;
        nbf[l + 3] = n + 3;
      }

      findNbFacets (k,     nn, nbf, Nf, fctp[i]);
      findNbFacets (k + 1, nn, nbf, Nf, fctp[i]);
      findNbFacets (k + 2, nn, nbf, Nf, fctp[i]);
      findNbFacets (k + 3, nn, nbf, Nf, fctp[i]);
    }

    Nf *= 4;
  }

  *nf = Nf;

  return fctp;
}

struct Vertex **createArray (int refinement)
{
  /* Creates and initializes an array of vertex arrays for multiple refinement levels. */
  struct Vertex **vtx = malloc (refinement * sizeof (struct Vertex*));

  for (int i = 0, n = 12; i < refinement; i++, n = 4 * (n - 2) + 2)

    vtx[i] = malloc (n * sizeof (struct Vertex));

  return vtx;
}

struct Facet **createFacets (int refinement)
{
  /* Creates and initializes an array of facet arrays for multiple refinement levels. */
  struct Facet **fct = malloc (refinement * sizeof (struct Facet*));

  for (int i = 0, n = 20; i < refinement; i++, n *= 4)

    fct[i] = malloc (n * sizeof (struct Facet));

  return fct;
}

void destroyFacetPoints (int refinement, struct FacetPoints **fctp)
{
  /* Frees memory allocated for the facet array. */
  for (int i = 0; i < refinement; i++)

    free (fctp[i]);

  free (fctp);
}

void destroyFacet (int refinement, struct Facet **fct)
{
  /* Frees memory allocated for the facet array. */
  for (int i = 0; i < refinement; i++)

    free (fct[i]);

  free (fct);
}

void destroyVertex (int refinement, struct Vertex **vtx)
{
  /* Frees memory allocated for the vertex array. */
  for (int i = 0; i < refinement; i++)

    free (vtx[i]);

  free (vtx);
}

void findNbFacets (int i, int nn, int nbf[nn], int nf,
                   struct FacetPoints fctp[nf])
{
  /* Identifies neighboring facets for a given facet based on dot product comparisons. */
  int n = 0;

  /* Iterate through potential neighboring facets */
  for (int ji = 0; ji < nn; ji++)
  {
    int j = nbf[ji];

    if (n == N_MAX_NGB_FACETS) break;

    /* Check if p1 of current facet matches any vertex of neighbor */
    if (dot (&fctp[i].p1, &fctp[j].p1) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p1, &fctp[j].p2) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p1, &fctp[j].p3) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    /* Check if p2 of current facet matches any vertex of neighbor */
    if (dot (&fctp[i].p2, &fctp[j].p1) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p2, &fctp[j].p2) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p2, &fctp[j].p3) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    /* Check if p3 of current facet matches any vertex of neighbor */
    if (dot (&fctp[i].p3, &fctp[j].p1) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p3, &fctp[j].p2) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }

    if (dot (&fctp[i].p3, &fctp[j].p3) > MIN_DOT_NODES)
    {
      fctp[i].nbf[n++] = j; continue;
    }
  }

  /* Store the total number of neighboring facets found */
  fctp[i].nnf = n;
}

void findNodes (int refinement,
                struct FacetPoints **fctp,
                struct Vertex **vtx,
                struct Facet **fct)
{
  /* Identifies and assigns node indices and coordinates based on facet connectivity. */
  for (int r = 0, nv = 12, nf = 20; r < refinement; r++)
  {
    bool (*NotInGroup)[N_CORNERS] = malloc (sizeof (bool[nf][N_CORNERS]));

    /* Initialize flags for each corner of every facet */
    for (int i = 0; i < nf; i++)
    {
      NotInGroup[i][0] = true;
      NotInGroup[i][1] = true;
      NotInGroup[i][2] = true;
    }

    struct IndexGroup (*Ig)[N_CORNERS] = malloc (sizeof (struct IndexGroup[nf][N_CORNERS]));

    int c = 0;

    /* Process each facet and its corners to assign nodes */
    for (int i = 0; i < nf; i++)

      for (int j = 0; j < N_CORNERS; j++)
      {
        fct[r][i].nnf = fctp[r][i].nnf;

        for (int k = 0; k < fctp[r][i].nnf; k++)

          fct[r][i].nbf[k] = fctp[r][i].nbf[k];

        fct[r][i].n_sph   = fctp[r][i].n_sph;
        fct[r][i].min_dot = fctp[r][i].min_dot;

        int n = 0;

        if (NotInGroup[i][j])
        {
          struct Triplet p;

          /* Select the appropriate vertex based on corner index */
          switch (j)
          {
            case 0:
              p = fctp[r][i].p1;
              break;

            case 1:
              p = fctp[r][i].p2;
              break;

            case 2:
              p = fctp[r][i].p3;
              break;
          }

          /* Find neighboring facets sharing this vertex */
          for (int k = 0; k < fctp[r][i].nnf; k++)
          {
            int index = fctp[r][i].nbf[k];

            if (n == N_MAX_NGB_TRIAG) break;

            if (dot (&p, &fctp[r][index].p1) > MIN_DOT_NODES)
            {
              Ig[i][j].id[n]   = 0;
              Ig[i][j].ii[n++] = index;

              NotInGroup[index][0] = false;

              continue;
            }

            if (dot (&p, &fctp[r][index].p2) > MIN_DOT_NODES)
            {
              Ig[i][j].id[n]   = 1;
              Ig[i][j].ii[n++] = index;

              NotInGroup[index][1] = false;

              continue;
            }

            if (dot (&p, &fctp[r][index].p3) > MIN_DOT_NODES)
            {
              Ig[i][j].id[n]   = 2;
              Ig[i][j].ii[n++] = index;

              NotInGroup[index][2] = false;

              continue;
            }
          }
        }

        if (n > 0)
        {
          struct Triplet p; p.x = 0; p.y = 0; p.z = 0;

          vtx[r][c].p.x = 0;
          vtx[r][c].p.y = 0;
          vtx[r][c].p.z = 0;

          /* Average coordinates from neighboring vertices */
          for (int k = 0; k < n; k++)
          {
            int ii = Ig[i][j].ii[k];
            int id = Ig[i][j].id[k];

            switch (id)
            {
              case 0:
                p = fctp[r][ii].p1;
                break;

              case 1:
                p = fctp[r][ii].p2;
                break;

              case 2:
                p = fctp[r][ii].p3;
                break;
            }

            vtx[r][c].p.x += p.x;
            vtx[r][c].p.y += p.y;
            vtx[r][c].p.z += p.z;
          }

          vtx[r][c].p.x /= n;
          vtx[r][c].p.y /= n;
          vtx[r][c].p.z /= n;

          if (fabs (vtx[r][c].p.x) < EPSILON) vtx[r][c].p.x = 0.0;
          if (fabs (vtx[r][c].p.y) < EPSILON) vtx[r][c].p.y = 0.0;
          if (fabs (vtx[r][c].p.z) < EPSILON) vtx[r][c].p.z = 0.0;

          /* Update facet indices with the new node */
          for (int k = 0; k < n; k++)
          {
            int ii = Ig[i][j].ii[k];
            int id = Ig[i][j].id[k];

            switch (id)
            {
              case 0:
                fct[r][ii].i1 = c;
                break;

              case 1:
                fct[r][ii].i2 = c;
                break;

              case 2:
                fct[r][ii].i3 = c;
                break;
            }
          }

          c++;
        }
      }

    /* Determine specific neighboring facet relationships */
    for (int i = 0; i < nf; i++)
    {
      int n = fct[r][i].nnf;

      for (int j = 0; j < n; j++)
      {
        int k = fct[r][i].nbf[j];

        if (k == i) continue;

        int c1 = 0, c2 = 0, c3 = 0;

        if ((fct[r][i].i1 == fct[r][k].i1) ||
            (fct[r][i].i1 == fct[r][k].i2) ||
            (fct[r][i].i1 == fct[r][k].i3)) c1++;

        if ((fct[r][i].i2 == fct[r][k].i1) ||
            (fct[r][i].i2 == fct[r][k].i2) ||
            (fct[r][i].i2 == fct[r][k].i3)) c2++;

        if ((fct[r][i].i3 == fct[r][k].i1) ||
            (fct[r][i].i3 == fct[r][k].i2) ||
            (fct[r][i].i3 == fct[r][k].i3)) c3++;

        if (c1 + c2 + c3 < 2) continue;

        if      (c1 == 1 && c2 == 1) fct[r][i].nbf_12 = k;
        else if (c1 == 1 && c3 == 1) fct[r][i].nbf_13 = k;
        else if (c2 == 1 && c3 == 1) fct[r][i].nbf_23 = k;
      }
    }

    nv = 4 * (nv - 2) + 2;
    nf = 4 * nf;

    free (NotInGroup);
    free (Ig);
  }
}

void findNbVertices (int refinement,
                     struct Facet **fct,
                     struct Vertex **vtx)
{
  /* Determines neighboring vertices for each vertex across refinement levels. */
  int nf = 20;
  int nv = 12;

  for (int r = 0; r < refinement; r++)
  {
    /* Process each facet to find neighboring vertices */
    for (int i = 0; i < nf; i++)
    {
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;

      int i1 = fct[r][i].i1;
      int i2 = fct[r][i].i2;
      int i3 = fct[r][i].i3;

      /* Check neighbors of each vertex of the current facet */
      for (int ji = 0; ji < fct[r][i].nnf; ji++)
      {
        int j = fct[r][i].nbf[ji];

        int j1 = fct[r][j].i1;
        int j2 = fct[r][j].i2;
        int j3 = fct[r][j].i3;

        if (i1 == j1)
        {
          if (notInList (j2, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j2;

          if (notInList (j3, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j3;
        }

        if (i1 == j2)
        {
          if (notInList (j1, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j1;

          if (notInList (j3, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j3;
        }

        if (i1 == j3)
        {
          if (notInList (j1, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j1;

          if (notInList (j2, n1, vtx[r][i1].nbv))

            vtx[r][i1].nbv[n1++] = j2;
        }

        if (i2 == j1)
        {
          if (notInList (j2, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j2;

          if (notInList (j3, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j3;
        }

        if (i2 == j2)
        {
          if (notInList (j1, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j1;

          if (notInList (j3, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j3;
        }

        if (i2 == j3)
        {
          if (notInList (j1, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j1;

          if (notInList (j2, n2, vtx[r][i2].nbv))

            vtx[r][i2].nbv[n2++] = j2;
        }

        if (i3 == j1)
        {
          if (notInList (j2, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j2;

          if (notInList (j3, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j3;
        }

        if (i3 == j2)
        {
          if (notInList (j1, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j1;

          if (notInList (j3, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j3;
        }

        if (i3 == j3)
        {
          if (notInList (j1, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j1;

          if (notInList (j2, n3, vtx[r][i3].nbv))

            vtx[r][i3].nbv[n3++] = j2;
        }
      }

      vtx[r][i1].nnv = n1;
      vtx[r][i2].nnv = n2;
      vtx[r][i3].nnv = n3;
    }

    nf *= 4; nv = 4 * (nv - 2) + 2;
  }
}

int findNextFacet (int i,
                   struct Triplet *p,
                   int nf, struct Facet fct[nf],
                   int nv, struct Vertex vtx[nv])
{
  /* Finds the next facet containing a projected point. */
  if (isInsideFacet (i, p, nf, fct, nv, vtx))

    return i;

  else

    /* Iterate through neighboring facets */
    for (int j = 0; j < fct[i].nnf; j++)
    {
      int jj = fct[i].nbf[j];

      if (jj == i) continue;

      if (isInsideFacet (jj, p, nf, fct, nv, vtx))

        return jj;
    }

  return -1;
}

int _findFacet (int refinement,
                int r, int i, int nf, int nv,
                struct Triplet p,
                struct Facet **fct,
                struct Vertex **vtx)
{
  /* Recursively searches for a facet containing a point at refinement level r. */
  if (dot (&fct[r][i].n_sph, &p) >= fct[r][i].min_dot)
  {
    struct Triplet p1 = vtx[r][fct[r][i].i1].p;
    struct Triplet p2 = vtx[r][fct[r][i].i2].p;
    struct Triplet p3 = vtx[r][fct[r][i].i3].p;

    struct Triplet q = project (&p, &p1, &p2, &p3);

    /* Baricentric coordinates */
    double u, v, w; baricentric (&q, &p1, &p2, &p3,
                                 &u, &v, &w);

    /* Check if point is inside facet */
    if (u + v + w > LEEWAY)

      return -1;

    else

      if (r + 1 == refinement)

        return i;

      else

        /* Recursively search for facet at next refinement level */
        for (int s = 0; s < 4; s++)
        {
          int index = _findFacet (refinement,
                                  r + 1,
                                  4 * i + s,
                                  4 * nf,
                                  4 * (nv - 2) + 2,
                                  p, fct, vtx);

          if (index >= 0) return index;
        }
  }

  return -1;
}

int findFacet (int refinement, 
               struct Triplet p,
               struct Facet **fct,
               struct Vertex **vtx)
{
  /* Finds the facet containing a projected point across refinement levels. */
  int nv = 12, nf = 20;

  for (int i = 0; i < nf; i++)
  {
    int index = _findFacet (refinement, 0, i, nf, nv, p, fct, vtx);

    if (index >= 0) return index;
  }

  return -1;
}

double interpolateElement (int i, struct Triplet *p,
                           double r1, double r2, double r3,
                           double r4, double r5, double r6,
                           double m1, double m2, double m3,
                           double m4, double m5, double m6,
                           int nf, struct Facet fct[nf],
                           int nv, struct Vertex vtx[nv])
{
  /* Interpolates a value within a prismatic element using weights. */
  double w1, w2, w3, w4, w5, w6;

  weights (p,
           r1, r2, r3, r4, r5, r6,
           &fct[i], nv, vtx,
           &w1, &w2, &w3, &w4, &w5, &w6);

  return m1 * w1 + m2 * w2 + m3 * w3
       + m4 * w4 + m5 * w5 + m6 * w6;
}
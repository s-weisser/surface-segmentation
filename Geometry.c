/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         Geometry.c
* AUTORS:       Sergej Rjasanow and Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Geometry.h"

#define LARGE_NUMBER 1e+16
#define BUFSIZE 4096

void cblas_dcopy(const long N, const double *X, const long incX, 
                 double *Y, const long incY);
void cblas_daxpy(const long N, const double alpha, const double *X,
                 const long incX, double *Y, const long incY);
double cblas_dnrm2(const long N, const double *X, const long incX);
void cblas_dscal(const long N, const double alpha, double *X, const long incX);

static void FindNodeNeighbours(long NumNodes, Node *Nodes, 
                               long NumElements, Triangle *Triangles) {
  long i,j,NumNbAll;
/*
  Count the neighbours 
*/
  for(i=0; i < NumNodes; i++)
    Nodes[i].NumNb=0;

  for(i=0; i < NumElements; i++)
    {
      Nodes[Triangles[i].Nodes[0]->GlNum].NumNb++;
      Nodes[Triangles[i].Nodes[1]->GlNum].NumNb++;
      Nodes[Triangles[i].Nodes[2]->GlNum].NumNb++;
    }
/*
  Number of neighbours
*/
  NumNbAll=0;
  for(i=0; i < NumNodes; i++)
    NumNbAll+=Nodes[i].NumNb;
/*
  Allocation of the memory for the neighbours for each node
*/
  Nodes[0].Nb=malloc(NumNbAll*sizeof(Triangle*));
  assert(Nodes[0].Nb != NULL);
/*
  Allocation of the memory for the neighbours for each node
*/
  Nodes[0].KNb=malloc(NumNbAll*sizeof(long));
  assert(Nodes[0].KNb != NULL);
/*
  Fill the field of neighbours
*/
  for(i=1; i < NumNodes; i++)
    {
      Nodes[i].Nb=Nodes[i-1].Nb+Nodes[i-1].NumNb;
      Nodes[i].KNb=Nodes[i-1].KNb+Nodes[i-1].NumNb;
      Nodes[i-1].NumNb=0;
    }
  Nodes[NumNodes-1].NumNb=0;
 
  for(i=0; i < NumElements; i++)
    { 
      j=Triangles[i].Nodes[0]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=0;

      j=Triangles[i].Nodes[1]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=1;

      j=Triangles[i].Nodes[2]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=2;
    }
}

static void FindTwo(long na,long nb,long *lista,long *listb,long *i1,long *i2) {
  long ia,ib,count;

  count=0;

  for (ia=0; ia < na; ia++)
    {
      for (ib=0; ib < nb; ib++)
	{
	  if (lista[ia] == listb[ib])
	    {
	      if (count == 0)
		{
		  *i1=lista[ia];
		  count=1;
		}
	      else if (count == 1)
		{
		  *i2=lista[ia];
		  count=2;
		}
	    }
	}

      if(count ==2) break;
    }
}

static double max(double a, double b) {
    if (a >= b) 
	return a; 
    else 
	return b;
}

static long FindElementNeighbours(long NumElements, Triangle *Triangles) {
  long i,i1=0,i2=0,ia,ib,nta,ntb;
  long *lista,*listb;
  Node *Na,*Nb;
/*
  All elements 
*/
  for(i=0; i < NumElements; i++)
    {
/*
  First edge
*/
      Na=Triangles[i].Nodes[0];
      Nb=Triangles[i].Nodes[1];

      nta=Na->NumNb;
      ntb=Nb->NumNb;

      lista=(long *)malloc(nta*sizeof(long));
      assert(lista != NULL);

      listb=(long *)malloc(ntb*sizeof(long));
      assert(listb != NULL);

      for (ia=0; ia < nta; ia++)
	lista[ia]=Na->Nb[ia]->GlNum;

      for (ib=0; ib < ntb; ib++)
	listb[ib]=Nb->Nb[ib]->GlNum;

      FindTwo(nta,ntb,lista,listb,&i1,&i2);

      free(lista); free(listb);

      if (i1 == i)
	Triangles[i].Nb[0]=Triangles+i2;
      else if (i2 == i)
	Triangles[i].Nb[0]=Triangles+i1;
      else
	return 1;
/*
  Second edge
*/
      Na=Triangles[i].Nodes[1];
      Nb=Triangles[i].Nodes[2];

      nta=Na->NumNb;
      ntb=Nb->NumNb;

      lista=(long *)malloc(nta*sizeof(long));
      assert(lista != NULL);

      listb=(long *)malloc(ntb*sizeof(long));
      assert(listb != NULL);

      for (ia=0; ia < nta; ia++)
	lista[ia]=Na->Nb[ia]->GlNum;

      for (ib=0; ib < ntb; ib++)
	listb[ib]=Nb->Nb[ib]->GlNum;

      FindTwo(nta,ntb,lista,listb,&i1,&i2);

      free(lista); free(listb);

      if (i1 == i)
	Triangles[i].Nb[1]=Triangles+i2;
      else if (i2 == i)
	Triangles[i].Nb[1]=Triangles+i1;
      else
	return 1;
/*
  Third edge
*/
      Na=Triangles[i].Nodes[2];
      Nb=Triangles[i].Nodes[0];

      nta=Na->NumNb;
      ntb=Nb->NumNb;

      lista=(long *)malloc(nta*sizeof(long));
      assert(lista != NULL);

      listb=(long *)malloc(ntb*sizeof(long));
      assert(listb != NULL);

      for (ia=0; ia < nta; ia++)
	lista[ia]=Na->Nb[ia]->GlNum;

      for (ib=0; ib < ntb; ib++)
	listb[ib]=Nb->Nb[ib]->GlNum;

      FindTwo(nta,ntb,lista,listb,&i1,&i2);

      free(lista); free(listb);

      if (i1 == i)
	Triangles[i].Nb[2]=Triangles+i2;
      else if (i2 == i)
	Triangles[i].Nb[2]=Triangles+i1;
      else
	return 1;
    }

  return 0;
}

static Triangle *FindOtherTriaNb(Node *node1, Node *node2, Triangle *tria, long *indexEdge) {
  long i;

  for (i=0; i<node1->NumNb; i++) 
    if ( node1->Nb[i] != tria &&
         (node2 == node1->Nb[i]->Nodes[0] || node2 == node1->Nb[i]->Nodes[1] || node2 == node1->Nb[i]->Nodes[2])) {
      *indexEdge = (node1->KNb[i]+2)%3; // local index of node2 in neighbouring triangle
      return node1->Nb[i];
    }

  return NULL;
}

static void FindEdgesAndNeighbours(long numNodes, Node *nodes, \
                                   long numEdges, Edge *edges, \
                                   long numTria, Triangle *triangles) {
  long i,j;
  long numEdgesFilled = 0;
  long indexEdge;
  Triangle *tria;

  /*
     Find edges and fill content.
  */

  for (i=0; i<numTria; ++i) {
    triangles[i].Edges[0] = NULL;
    triangles[i].Edges[1] = NULL;
    triangles[i].Edges[2] = NULL;
  }

  // go over all triangles and add edges if not already done
  for (i=0; i<numTria; i++) {
    for (j=0; j<3; ++j) {
      if (triangles[i].Edges[j] == NULL) {
        tria = FindOtherTriaNb(triangles[i].Nodes[j], triangles[i].Nodes[(j+1)%3], &(triangles[i]), &indexEdge);
        assert(tria != NULL);

        edges[numEdgesFilled].GlNum = numEdgesFilled;
        edges[numEdgesFilled].Nodes[0] = triangles[i].Nodes[j];
        edges[numEdgesFilled].Nodes[1] = triangles[i].Nodes[(j+1)%3];
        edges[numEdgesFilled].NumNb = 2;
        edges[numEdgesFilled].Nb[0] = &(triangles[i]);
        edges[numEdgesFilled].Nb[1] = tria;

        triangles[i].Edges[j]  = &(edges[numEdgesFilled]);
        tria->Edges[indexEdge] = &(edges[numEdgesFilled]);

        numEdgesFilled++;
      }
    }
  }
}

static void dvecpr(double X[3],double Y[3],double Z[3])
{
  Z[0]=X[1]*Y[2]-X[2]*Y[1];
  Z[1]=X[2]*Y[0]-X[0]*Y[2];
  Z[2]=X[0]*Y[1]-X[1]*Y[0];
}

static long ElementInfo(long NumElements, Triangle *Triangles) {
  double U12[3],U23[3],U31[3],Q,d12,d23,d31,d;
  long i;
  double *X1,*X2,*X3;
  Triangle *TR;

  for(i=0; i < NumElements; i++)
    {
      TR=Triangles+i;

      X1=TR->Nodes[0]->X;
      X2=TR->Nodes[1]->X;
      X3=TR->Nodes[2]->X;

      cblas_dcopy(3,X2,1,U12,1);
      cblas_daxpy(3,-1.0,X1,1,U12,1);
      d12=cblas_dnrm2(3,U12,1);

      cblas_dcopy(3,X3,1,U23,1);
      cblas_daxpy(3,-1.0,X2,1,U23,1);
      d23=cblas_dnrm2(3,U23,1);

      cblas_dcopy(3,X1,1,U31,1);
      cblas_daxpy(3,-1.0,X3,1,U31,1);
      d31=cblas_dnrm2(3,U31,1);

      d=max(max(d12,d23),d31);

      dvecpr(U31,U12,TR->Norm);
      Q=cblas_dnrm2(3,TR->Norm,1);
      if (d*d >= Q*LARGE_NUMBER)
	  return i+1;
      else
	  cblas_dscal(3,1.0/Q,TR->Norm,1);
    }
  return 0;
}

int ReadMesh(char *filename, Mesh *mesh) {
  long i, j;
  char buf[BUFSIZE];  
  char *ec;
  int  res;
  FILE *GeometryFile;
  long GN, numTags;
  long idx_node0, idx_node1, idx_node2, idx_triangle;

  long numNodes, numEdges, numElements, numTriangles;
  Node *nodes;
  Edge *edges;
  Triangle *triangles;

/* Open the file with the geometry */
  GeometryFile = fopen(filename, "r"); 
  if ( GeometryFile == NULL ) {
    fprintf(stderr, "Could not open %s\n", filename);
    return 1;
  }

/* Check file format */
  ec = fgets(buf, BUFSIZE, GeometryFile);
  if (ec==NULL)
    return 4;

  if ( !strstr(buf, "$MeshFormat") ) {
    fprintf(stderr, "Not suported mesh format in %s.\nUse .msh from gmsh in ASCII format.\n", filename);
    return 2;
  }

  ec = fgets(buf, BUFSIZE, GeometryFile);
  if (ec==NULL)
    return 4;

  if ( !strstr(buf, "2.2 0 8") ) {
    fprintf(stderr, "Not suported mesh format in %s.\nUse .msh from gmsh in ASCII format.\n", filename);
    return 2;
  }

/* Search for the beginning of the node data */
  do {
    ec = fgets(buf, BUFSIZE, GeometryFile);
    if (ec==NULL)
      return 4;
  } while(strstr(buf, "$Nodes") == 0);
  ec = fgets(buf, BUFSIZE, GeometryFile);
  if (ec==NULL)
    return 4;

  sscanf(buf, "%ld", &numNodes);
  nodes = calloc(numNodes,sizeof(Node));
  if (nodes == NULL) 
    return 3;

/* Read nodes */
  for (i=0; i<numNodes; i++) {
    ec = fgets(buf, BUFSIZE, GeometryFile);
    if (ec==NULL)
      return 4;

    sscanf(buf, "%ld %lf %lf %lf", &GN, &nodes[i].X[0], &nodes[i].X[1], &nodes[i].X[2]);
    nodes[i].GlNum = GN-1;
  }

/* Search for the beginning of the element data */
  do {
    ec = fgets(buf, BUFSIZE, GeometryFile);
    if (ec==NULL)
      return 4;
  } while(strstr(buf, "$Elements") == 0);
  ec = fgets(buf, BUFSIZE, GeometryFile);
  if (ec==NULL)
    return 4;

  sscanf(buf, "%ld", &numElements);

/* Count number of triangles in surface mesh */
  numTriangles = 0;
  for(i=0; i<numElements; i++) {
    ec = fgets(buf, BUFSIZE, GeometryFile);
    if (ec==NULL)
      return 4;

    if ( sscanf(buf, "%ld 2 %ld", &GN, &numTags) == 2 )
      numTriangles++;
  }

/* Search for the beginning of the element data */
  rewind(GeometryFile);
  do {
    ec = fgets(buf, BUFSIZE, GeometryFile);
    if (ec==NULL)
      return 4;
  } while(strstr(buf, "$Elements") == 0);
  ec = fgets(buf, BUFSIZE, GeometryFile);
  if (ec==NULL)
    return 4;

  triangles = calloc(numTriangles,sizeof(Triangle));
  if (triangles == NULL) 
    return 3;

/* Read elements, i.e. triangles */
  idx_triangle = 0;
  for(i=0; i<numElements; i++) {
    if ( fscanf(GeometryFile, "%*d 2 %ld", &numTags) == 1 ) { // element is a triangle?
      for (j=0; j<numTags; j++) 
        res = fscanf(GeometryFile, "%*d");

      res = fscanf(GeometryFile, "%ld %ld %ld", &idx_node0, &idx_node1, &idx_node2);
      if (res!=3)
        return 4;

      triangles[idx_triangle].GlNum = idx_triangle;
      triangles[idx_triangle].Nodes[0] = &nodes[idx_node0-1];
      triangles[idx_triangle].Nodes[1] = &nodes[idx_node1-1];
      triangles[idx_triangle].Nodes[2] = &nodes[idx_node2-1];
      idx_triangle++;
    }
    else {
      ec = fgets(buf, BUFSIZE, GeometryFile);
      if (ec==NULL)
        return 4;
    }
  }

  fclose(GeometryFile);

  numEdges = 1.5*numTriangles;
  edges = calloc(numEdges,sizeof(Edge));
  if (edges == NULL)
    return 3;

  mesh->NumNodes = numNodes;
  mesh->NumEdges = numEdges;
  mesh->NumElements = numTriangles;
  mesh->Nodes = nodes;
  mesh->Edges = edges;
  mesh->Triangles = triangles;

  return 0;
}

void FreeGeometry(Mesh *mesh) {
  free(mesh->Nodes[0].Nb);
  free(mesh->Nodes[0].KNb);
  free(mesh->Nodes);
  free(mesh->Edges);
  free(mesh->Triangles);
  mesh->Nodes = NULL;
  mesh->Edges = NULL;
  mesh->Triangles = NULL;
}

int SetupGeometry(Mesh mesh) {
  int err = 0;
  FindNodeNeighbours(mesh.NumNodes, mesh.Nodes, mesh.NumElements, mesh.Triangles);
  FindElementNeighbours(mesh.NumElements, mesh.Triangles);
  err = ElementInfo(mesh.NumElements, mesh.Triangles);
  FindEdgesAndNeighbours(mesh.NumNodes, mesh.Nodes, 
                         mesh.NumEdges, mesh.Edges, 
                         mesh.NumElements, mesh.Triangles);
  return err;
}

int WriteMesh_msh(char *filename, Mesh mesh) {
  long i;
  FILE *file;

  file = fopen(filename, "w"); 
  if ( file == NULL ) {
    fprintf(stderr, "Could not open %s\n", filename);
    return 1;
  }

  fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

  /* write nodes */
  fprintf(file, "$Nodes\n%ld\n", mesh.NumNodes);
  for (i=0; i<mesh.NumNodes; i++) {
    fprintf(file, "%ld %e %e %e\n", i+1, mesh.Nodes[i].X[0], mesh.Nodes[i].X[1], mesh.Nodes[i].X[2]);
  }
  fprintf(file, "$EndNodes\n");

  /* write elements, i.e. triangles */
  fprintf(file, "$Elements\n%ld\n", mesh.NumElements);
  for (i=0; i<mesh.NumElements; i++) {
    fprintf(file, "%ld 2 0 %ld %ld %ld\n", i+1, mesh.Triangles[i].Nodes[0]->GlNum+1, 
                                                mesh.Triangles[i].Nodes[1]->GlNum+1, 
                                                mesh.Triangles[i].Nodes[2]->GlNum+1);
  }
  fprintf(file, "$EndElements\n");

  /* write element data */
  fprintf(file, "$ElementData\n");
  fprintf(file, "1\n\"segment number\"\n1\n0.0\n3\n0\n1\n%ld\n", mesh.NumElements);
  for (i=0; i<mesh.NumElements; i++) {
    fprintf(file, "%ld %ld\n", i+1, mesh.Triangles[i].SegmentNum);
  }
  fprintf(file, "$EndElementData\n");

  fclose(file);

  return 0;
}

int WriteMesh_vtk(char *filename, Mesh mesh) {
  long i;
  FILE *file;

  file = fopen(filename, "w"); 
  if ( file == NULL ) {
    fprintf(stderr, "Could not open %s\n", filename);
    return 1;
  }

  /* write header */
  fprintf(file, "# vtk DataFile Version 3.0\nSurface Segmentation by S. Weisser\nASCII\n");

  fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

  /* write nodes */
  fprintf(file, "POINTS %ld double\n", mesh.NumNodes);
  for (i=0; i<mesh.NumNodes; i++) {
    fprintf(file, "%e %e %e\n", mesh.Nodes[i].X[0], mesh.Nodes[i].X[1], mesh.Nodes[i].X[2]);
  }

  /* write elements */
  fprintf(file, "CELLS %ld %ld\n", mesh.NumElements, 4 * mesh.NumElements);
  for (i=0; i<mesh.NumElements; i++) {
    fprintf(file, "3 %ld %ld %ld\n", mesh.Triangles[i].Nodes[0]->GlNum, 
                                     mesh.Triangles[i].Nodes[1]->GlNum, 
                                     mesh.Triangles[i].Nodes[2]->GlNum);
  }

  /* specify elements as triangles */
  fprintf(file, "CELL_TYPES %ld\n", mesh.NumElements);
  for (i=0; i<mesh.NumElements; i++) {
    fprintf(file, "5\n");
  }

  /* write element data */
  fprintf(file, "CELL_DATA %ld\n", mesh.NumElements);
  fprintf(file, "SCALARS segment_number int\nLOOKUP_TABLE default\n");
  for (i=0; i<mesh.NumElements; i++) {
    fprintf(file, "%ld\n", mesh.Triangles[i].SegmentNum);
  }

  fclose(file);
  
  return 0;
}

#undef LARGE_NUMBER
#undef BUFSIZE

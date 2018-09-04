/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         Geometry.h
* AUTORS:       Sergej Rjasanow and Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

/**************** declaration of variables and structs ***********************/

typedef struct sNode {
  long GlNum;
  double X[3];
  long NumNb;
  struct sTriangle **Nb;
  long  *KNb;
} Node;

typedef struct sEdge {
  long GlNum;
  Node *Nodes[2];
  long NumNb;
  struct sTriangle *Nb[2];
} Edge;

typedef struct sTriangle {
  long GlNum;
  long SegmentNum;
  Node *Nodes[3];
  Edge *Edges[3];
  double Norm[3];
  struct sTriangle *Nb[3];
} Triangle;

typedef struct sMesh {
  long NumNodes;
  Node *Nodes;
  long NumEdges;
  Edge *Edges;
  long NumElements;
  Triangle *Triangles;
} Mesh;

/**************** declaration of global functions ***************************/
int ReadMesh(char *filename, Mesh *mesh);
int SetupGeometry(Mesh mesh);
void FreeGeometry(Mesh *mesh);
int WriteMesh_msh(char *filename, Mesh mesh);
int WriteMesh_vtk(char *filename, Mesh mesh);

#endif

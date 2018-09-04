/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         PartitionBoundaryMesh.h
* AUTORS:       Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#ifndef PARTITIONBOUNDARYMESH_H
#define PARTITIONBOUNDARYMESH_H

#include "Geometry.h"

/*
 This functions generates a partitioning of the boundary mesh
 which is given in the parameter 'mesh'. For this reason geometric
 edges are first detected and then the partitions are formed such
 that they do not contain these edges. The partitioning follows
 the surface segmentation discussed in

 S. Rjasanow, S. Weißer: ACA improvement by surface segmentation, 2018 (submitted)

 Edges in the mesh are detected as geometric edges if the dihedral 
 angle of the neighbouring triangles exceeds a threshold angle.
 The input parameter of the functions relates to this angle by

    threshold = cos( angle )

 The triangles in the mesh are assigned to segments by setting
 the value 'mesh.Triangles[i].SegmentNum'. 
 The function returns the number of generated segments.
*/
long PartitionBoundaryMesh(double threshold, Mesh mesh);

#endif

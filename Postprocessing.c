/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         Postprocessing.c
* AUTORS:       Steffen Weisser
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
#include "Postprocessing.h"
#include "Geometry.h"

double cblas_ddot(const long N, const double *X, const long incX,
                  const double *Y, const long incY);

long UnifySegments(long numSegments, double threshold, Mesh mesh) {
  long i, j, id;
  short *matrix, unified;
  long num1, num2;
  long matrixDim;
  long *newSegmentIDs;
  double cos_angle;
  Edge *edge;

  matrixDim = numSegments;
  matrix = calloc(matrixDim*matrixDim,sizeof(short));
  newSegmentIDs = calloc(numSegments,sizeof(double));
  assert(newSegmentIDs!=NULL);

  for (i=0; i<numSegments; i++)
    newSegmentIDs[i] = i;

  // setup matrix indicating admissible segments for unification
  for (i=0; i<mesh.NumEdges; i++) {
    edge = &(mesh.Edges[i]);

    if (edge->Nb[0]->SegmentNum != edge->Nb[1]->SegmentNum) {
      num1 = edge->Nb[0]->SegmentNum;
      num2 = edge->Nb[1]->SegmentNum;

      cos_angle = cblas_ddot(3, edge->Nb[0]->Norm, 1, edge->Nb[1]->Norm, 1);

      if (cos_angle < threshold) { // geometric edge?
        matrix[num1*matrixDim+num2] = 1;
        matrix[num2*matrixDim+num1] = 1;
      }
      else {
        if (matrix[num1*matrixDim+num2]==0) {
          matrix[num1*matrixDim+num2] = -1;
          matrix[num2*matrixDim+num1] = -1;
        }
      } 
    }
  }

  do {
    unified = 0;
    for (num1=0; num1<matrixDim; num1++) {
      if (num1==newSegmentIDs[num1]) {
        for (num2=num1+1; num2<matrixDim; num2++) {
          if (num2==newSegmentIDs[num2]) {
            // two segments found which might be unified
            // check if they are admissible
            if (matrix[num1*matrixDim+num2]==-1) {
              // unify the segments with numbers num1 and num2 and update matrix
              for (i=0; i<matrixDim; i++) {
                if (matrix[num2*matrixDim+i]==1) {
                  matrix[num1*matrixDim+i] = 1;
                } else if (matrix[num2*matrixDim+i]==-1
                           && matrix[num1*matrixDim+i]!=1) {
                  matrix[num1*matrixDim+i] = -1;
                }
                // update newSegmentIDs
                if (newSegmentIDs[i]==num2)
                  newSegmentIDs[i] = num1;
              }
              unified++;
              numSegments--;

              break;
            }
          }
        }
      }
    }
  } while (unified);

  // adjust numbering from 0 to numSegments
  id = 0;
  for (i=0; i<matrixDim; i++) {
    if (newSegmentIDs[i]==i) {
      for (j=0; j<matrixDim; j++) {
        if (newSegmentIDs[j]==i)
          newSegmentIDs[j] = id;
      }
      id++;
    }
  }
  
  // update SegmentNum in triangles
  for (i=0; i<mesh.NumElements; i++)
    mesh.Triangles[i].SegmentNum = newSegmentIDs[mesh.Triangles[i].SegmentNum];

  free(matrix);
  free(newSegmentIDs);

  return numSegments;
}

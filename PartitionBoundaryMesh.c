/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         PartitionBoundaryMesh.c
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
#include "PartitionBoundaryMesh.h"
#include "Geometry.h"
#include "List.h"

#define LENGTHPATCHESARRAY 64

double cblas_ddot(const long N, const double *X, const long incX,
                  const double *Y, const long incY);

static void followEdge(Edge *edge, Node *node,
                       double threshold, short *edgeTouched, 
                       long *numSegments, List ***segments);

long PartitionBoundaryMesh(double threshold, Mesh mesh) {
  long i,j,k;
  double cos_angle;
  long numSegments = 0, numSegmentsFinished, numTriaInSegment;
  short *edgeTouched;
  Triangle *tria;

  // for each segment we realize a list of triangles
  List **segments;

  long numEdges       = mesh.NumEdges;
  long numTria        = mesh.NumElements;
  Edge *edges         = mesh.Edges;
  Triangle *triangles = mesh.Triangles;

  // initialize SegmentNum as not assigned
  for (i=0; i<numTria; i++)
    triangles[i].SegmentNum = -1;

  // allocate auxiliary data structures
  edgeTouched = (short*)calloc(numEdges,sizeof(short));
  segments = (List**)calloc(LENGTHPATCHESARRAY,sizeof(List*));
  assert(edgeTouched != NULL && segments != NULL);

  // degect geometric edges (Algorithm 2)
  for (i=0; i<numEdges; i++) {
    if (edgeTouched[i]) 
      continue;

    cos_angle = cblas_ddot(3, edges[i].Nb[0]->Norm, 1, edges[i].Nb[1]->Norm, 1);

    if (cos_angle < threshold) { // geometric edge?
      // mark edge as touched/processed
      edgeTouched[i] = 1;

      // assign neighbouring triangles a new SegmentNum
      if (numSegments+2 > LENGTHPATCHESARRAY) {
        segments = (List**)realloc(segments,(numSegments+2)*sizeof(List*));
        assert(segments != NULL);
      }

      // assign left triangle (new segment number)
      edges[i].Nb[0]->SegmentNum = numSegments;

      // initialize list for new segment
      segments[numSegments] = (List*)malloc(sizeof(List));
      assert(segments[numSegments]!=NULL);
      ListInit(segments[numSegments]);

      // add left triangle to the list
      ListAdd(segments[numSegments],(void*)edges[i].Nb[0]);
      numSegments++;

      // assign right triangle (new segment number)
      edges[i].Nb[1]->SegmentNum = numSegments;

      // initialize list for new segment
      segments[numSegments] = (List*)malloc(sizeof(List));
      assert(segments[numSegments]!=NULL);
      ListInit(segments[numSegments]);

      // add right triangle to the list
      ListAdd(segments[numSegments],(void*)edges[i].Nb[1]);
      numSegments++;

      // follow edge in both directions
      followEdge(&(edges[i]), edges[i].Nodes[0], threshold, edgeTouched, &numSegments, &segments);
      followEdge(&(edges[i]), edges[i].Nodes[1], threshold, edgeTouched, &numSegments, &segments);
    }
  }

  free(edgeTouched);

  if (numSegments > 0) {
    // complete the segments (Algorithm 4)
    do {
      numSegmentsFinished = 0;
      for (i=0; i<numSegments; i++) {
        numTriaInSegment = ListGetLength(segments[i]);
        if (numTriaInSegment == 0) {
          // only iterate over segments which still have
          // the possiblity to be enlarged
          numSegmentsFinished++;
          continue;
        }

        for (j=0; j<numTriaInSegment; j++) {
          tria = (Triangle*) ListGetFirst(segments[i]);

          // assign neighbouring triangles to the segment
          // if not already assigned to some segment
          for (k=0; k<3; k++) {
            if (tria->Nb[k]->SegmentNum < 0) {
              tria->Nb[k]->SegmentNum = i;
              ListAdd(segments[i],(void*)tria->Nb[k]);
            }
          }

          // delete processed triangles in segment list in order 
          // to iterate only over free boundary of segment
          ListDel(segments[i],(void*)tria);
        }
      }
    } while(numSegments != numSegmentsFinished);

    for (i=0; i<numSegments; i++) 
      free(segments[i]);
  }
  else {
    // there was no geometric edge and thus we have
    // only one segment containing the whole boundary
    for (i=0; i<numTria; i++)
      triangles[i].SegmentNum = 0;

    numSegments = 1;
  }

  free(segments);

  return numSegments;
}

/* 
 This function gets a node, an edge and a triangle, where the node is 
 an endpoint of the edge and the edge lies in the boundary of the triangle.
 The pointers edge and tria are changed in such a way that after the 
 function call edge points to the second edge of the triangle containing
 the node and tria points to the neighbouring triangle of the new edge.
*/
static void getOtherNeighbours(Node *node, Edge **edge, Triangle **tria) {
  int i;
  Edge *nextEdge;

  for (i=0; i<3; i++) {
    nextEdge = (*tria)->Edges[i];
    if (nextEdge != *edge && ( nextEdge->Nodes[0] == node || nextEdge->Nodes[1] == node))
      break;
  }

  if (i==3) {
    *edge = NULL;
    *tria = NULL;
    return;
  }

  if (nextEdge->Nb[0] == *tria)
    *tria = nextEdge->Nb[1];
  else
    *tria = nextEdge->Nb[0];

  *edge = nextEdge;
}

// follow edge to the direction of the node (Algorithm 3)
static void followEdge(Edge *edge, Node *node,
                       double threshold, short *edgeTouched, 
                       long *numSegments, List ***segments) {

  double cos_angle0, cos_angle1;
  Edge *nextEdge0, *nextEdge1, *lastEdge;
  Triangle *tria, *nextTria0, *nextTria1;
  List edgeList;

  // follow the geometric edge
  do {
    // initialize list to remember edges
    ListInit(&edgeList);

    // setup edgeList and initialize neighbouring triangles
    nextEdge0 = nextEdge1 = edge;
    nextTria0 = edge->Nb[0];
    nextTria1 = edge->Nb[1];
    cos_angle0 = cos_angle1 = 1;

    // process neighbourhood of node starting from edge
    // alternately in clockwise and anti clockwise direction
    // until a geometric edge is found each
    do {
      // check neighbourhood of node next to edge->Nb[0]
      if (cos_angle0 >= threshold) {
        tria = nextTria0;
        getOtherNeighbours(node, &nextEdge0, &nextTria0);
        
        cos_angle0 = cblas_ddot(3, tria->Norm, 1, nextTria0->Norm, 1);
  
        // in case that edge is not a geometric edge assign the 
        // triangle to the segment if it is not already assigned
        if (cos_angle0 >= threshold && nextTria0->SegmentNum < 0) {
          nextTria0->SegmentNum = tria->SegmentNum;
          ListAdd((*segments)[tria->SegmentNum],(void*)nextTria0);
        }
      }

      if (nextEdge0 == nextEdge1) 
        break;
      
      // check neighbourhood of node next to edge->Nb[1]
      if (cos_angle1 >= threshold) {
        tria = nextTria1;
        getOtherNeighbours(node, &nextEdge1, &nextTria1);
        
        cos_angle1 = cblas_ddot(3, tria->Norm, 1, nextTria1->Norm, 1);
  
        // in case that edge is not a geometric edge, assign the
        // triangle to the segment if it is not already assigned
        if (cos_angle1 >= threshold && nextTria1->SegmentNum < 0) {
          nextTria1->SegmentNum = tria->SegmentNum;
          ListAdd((*segments)[tria->SegmentNum],(void*)nextTria1);
        }
      }

      if (nextEdge0 == nextEdge1) 
        break;

    } while (cos_angle0 >= threshold && cos_angle1 >= threshold);

    if (nextEdge0 == nextEdge1 && cos_angle0 < threshold) { // circle closed and one edge to follow
      if (!edgeTouched[nextEdge0->GlNum]) {
        edgeTouched[nextEdge0->GlNum] = 1;
        ListAdd(&edgeList,(void*)nextEdge0);
      }
    } else if (nextEdge0 != nextEdge1) { // i.e., cos_angle0 < threshold && cos_angle1 < threshold
      // process remaining neighbouhood of node
      // close circle from direction edge->Nb[0]
      if (nextEdge0->Nb[0] == nextTria0)
        tria = nextEdge0->Nb[1];
      else
        tria = nextEdge0->Nb[0];

      do {
        if (!edgeTouched[nextEdge0->GlNum] && cos_angle0 < threshold) {
          // geometric edge detected, remember for later
          edgeTouched[nextEdge0->GlNum] = 1;
          ListAdd(&edgeList,(void*)nextEdge0);
        }

        if (nextTria0->SegmentNum < 0 && cos_angle0 < threshold) {
          // since there is a geometric edge, assign new segment number to the triangle
          nextTria0->SegmentNum = *numSegments;

          // initialize list for new segment
          if (*numSegments+1 > LENGTHPATCHESARRAY) {
            *segments = (List**)realloc(*segments,(*numSegments+1)*sizeof(List*));
            assert(*segments != NULL);
          }
          (*segments)[*numSegments] = (List*)malloc(sizeof(List));
          assert((*segments)[*numSegments]!=NULL);
          ListInit((*segments)[*numSegments]);
          ListAdd((*segments)[*numSegments],(void*)nextTria0);
          (*numSegments)++;
        } else if (nextTria0->SegmentNum < 0) {
          // since there is no geometric edge, assign previous segment number to the triangle
          nextTria0->SegmentNum = tria->SegmentNum;
          ListAdd((*segments)[tria->SegmentNum],(void*)nextTria0);
        }

        tria = nextTria0;
        getOtherNeighbours(node, &nextEdge0, &nextTria0);
        cos_angle0 = cblas_ddot(3, tria->Norm, 1, nextTria0->Norm, 1);
      } while (nextEdge0 != nextEdge1); 

      if (!edgeTouched[nextEdge0->GlNum] && cos_angle0 < threshold) {
        edgeTouched[nextEdge0->GlNum] = 1;
        ListAdd(&edgeList,(void*)nextEdge0);
      }
    }

    // follow each edge which was added to the edgeList
    lastEdge = (Edge*) ListGetLast(&edgeList);
    edge = (Edge*) ListGetFirst(&edgeList);

    if (edge != NULL) {
      if (edge->Nodes[0] == node)
        node = edge->Nodes[1];
      else
        node = edge->Nodes[0];

      // last Edge in endgeList will be treated with the outer loop
      while (edge != lastEdge) {
        // follow the first edges with recursion
        followEdge(edge, node, threshold, edgeTouched, numSegments, segments);

        edge = (Edge*) ListGetNext(&edgeList);

        if (edge->Nodes[0] == node)
          node = edge->Nodes[1];
        else
          node = edge->Nodes[0];
      }
  
      ListFree(&edgeList);
    }
  } while(edge!=NULL);

}


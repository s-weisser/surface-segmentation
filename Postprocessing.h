/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         Postprocessing.h
* AUTORS:       Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include "Geometry.h"

/*
 This functions reduces the number of segments by unifying 
 neighbouring segments that only share boundaries in the smooth 
 part of the surface. 

 Therefore the input parameters are the current number of segments,
 the threshold parameter in order detect geometric edges and the
 mesh. The values 'mesh.Triangles[i].SegmentNum' are updated.
 The return value of the function is the new number of segments. 
*/
long UnifySegments(long numSegments, double threshold, Mesh mesh);

#endif

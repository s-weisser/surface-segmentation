/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         SegmentBM.c
* AUTORS:       Sergej Rjasanow and Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "Geometry.h"
#include "PartitionBoundaryMesh.h"
#include "Postprocessing.h"

#define BUFSIZE 4096

void read_parameters(int, char*[],double*,char*,char*,short*);

int main(int argc, char *argv[]) {
  long NumSegments = 0;
  double threshold, angle;
  char prefix_output_filename[BUFSIZE];
  char filename[BUFSIZE];
  char *file_extension;
  short postprocessing;
  Mesh mesh;

  printf("\n");
  printf("Surface Segmentation\n");
  printf("by Sergej Rjasanow and Steffen Weisser\n");
  printf("Saarland University, Germany, 2018\n");
  printf("Copyright (C) 2018 Steffen Weisser - All Rights Reserved\n\n");

  read_parameters(argc,argv,&threshold,prefix_output_filename,filename,&postprocessing);

  file_extension = prefix_output_filename + strlen(prefix_output_filename) - 4;

  printf("Reading geometry from mesh file \"%s\"\n", filename);

  if ( ReadMesh(filename, &mesh) ) {
    fprintf(stderr, "ERROR while reading mesh.\n");
    return 1;
  }

  if ( SetupGeometry(mesh) ) {
    fprintf(stderr, "ERROR while generating geometric relations.\n");
    return 2;
  }

  printf("  number of nodes:\t%6ld\n  number of triangles:\t%6ld\n", mesh.NumNodes, mesh.NumElements);
  
  angle = acos(threshold);
  printf("\nPerforming surface segmentation:\n");
  printf("  threshold as angle:\t\t%5.2f (i.e. %5.2f degree)\n", angle, 180*M_1_PI*angle);

  // generate surface segmentation
  NumSegments = PartitionBoundaryMesh(threshold, mesh);
  printf("  detected surface segments:\t%5ld\n", NumSegments);

  if (postprocessing) {
    printf("\nPerforming postprocessing to reduce number of segments\n");
  
    NumSegments = UnifySegments(NumSegments, threshold, mesh);

    printf("  new number of segments:\t%5ld\n", NumSegments);
  }

  if ( strcmp(file_extension,".vtk") ) {
    if ( strcmp(file_extension,".msh") )
      sprintf(filename, "%s.msh", prefix_output_filename);
    else
      strncpy(filename, prefix_output_filename, BUFSIZE);

    printf("\nWriting output into \"%s\"\n", filename);
    if ( WriteMesh_msh(filename, mesh) ) {
      fprintf(stderr, "ERROR while writing msh mesh.\n");
      return 1;
    }
  }

  if ( strcmp(file_extension,".msh") ) {
    if ( strcmp(file_extension,".vtk") )
      sprintf(filename, "%s.vtk", prefix_output_filename);
    else
      strncpy(filename, prefix_output_filename, BUFSIZE);

    printf("\nWriting output into \"%s\"\n", filename);
    if ( WriteMesh_vtk(filename, mesh) ) {
      fprintf(stderr, "ERROR while writing msh mesh.\n");
      return 1;
    }
  }

  FreeGeometry(&mesh);
  
  puts("");

  return 0;
}

void read_parameters(int argc, char *argv[],
                     double *threshold, char *fileprefix, char *filename, short *postprocessing) {

  int c;
  double angle;
  const char *optstring = "a:A:t:o:ph";

  char *help = "Command line arguments:\n"
               "-a\t threshold parameter as angle between 0 and Pi/2 (default: Pi/3)\n"
               "-A\t threshold parameter as angle in degree between 0 and 90 (default: 60)\n"
               "-t\t threshold parameter as cos(angle) between 0 and 1 (default: 0.5)\n"
               "-o\t prefix of output files (default: out)\n\t use, i.e., out.vtk to generate only vtk file\n"
               "-p\t postprocessing to reduce number of segments\n"
               "-h\t show help message\n\n";

  *threshold = 0.5;
  strcpy(fileprefix, "out");
  *postprocessing = 0;
  
  while((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
      case 'a': 
        angle = strtod(optarg,NULL); 
        if ( angle>=0 && angle<=M_PI_2 )
          *threshold = cos(strtod(optarg,NULL)); 
        else 
          *threshold = -1;
        break;
      case 'A': 
        angle = strtod(optarg,NULL); 
        if ( angle>=0 && angle<=90 )
          *threshold = cos(M_PI_2/90. * strtod(optarg,NULL)); 
        else 
          *threshold = -1;
        break;
      case 't': *threshold = strtod(optarg,NULL); break;
      case 'o': strcpy(fileprefix, optarg); break;
      case 'p': *postprocessing = 1; break;
      case 'h': fputs(help,stdout); exit(0); break;
      default:  fputs(help,stderr); exit(-1); break;
    }
  }

  if ( *threshold<0 || *threshold>1 ) {
    fprintf(stderr, "ERROR: Threshold parameter out of range.\n\n");
    fputs(help,stderr); 
    exit(-1);
  }

  if (argv[optind]!=0) {
    strcpy(filename, argv[optind]);
  }
  else {
    fprintf(stderr, "USAGE: %s [atoph] meshfile.msh\n\n", argv[0]);
    fputs(help,stderr);
    exit(-1);
  }
}


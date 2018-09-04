#**************************************************************
#
# PROJECT:      Surface Segmentation
# FILE:         Makefile
# AUTORS:       Steffen Weisser
# AFFILIATION:  Saarland University, Germany
# DATE:         2018
#
# Copyright (C) 2018 Steffen Weisser - All Rights Reserved
# You may use, distribute and modify this code under the
# terms of the MIT license.
#
#**************************************************************

.SUFFIXES:	.c .o

##### Main executable ##################################################

MAIN 		= SegmentBM

##### Compilers ########################################################

CC              = gcc
LINKER          = gcc

##### Compiler and linker options ######################################

COP             = -Wall -O2 
LINKOP          = -lopenblas -lm 

##### source files #####################################################

SRC =		\
		SegmentBM.c \
		\
		Geometry.c \
                List.c \
                PartitionBoundaryMesh.c \
                Postprocessing.c \

OBJ = ${SRC:.c=.o}

########################################################################
##### Common part, please, do not edit ! ###############################
########################################################################

all:		$(MAIN) 

$(MAIN):	$(OBJ) Makefile
		$(LINKER) $(OBJ) -o $(MAIN) $(LINKOP)

.c.o:		
		$(CC) $(COP) -c $<

.PHONY:		clean, realclean, depend

clean:
	rm -f $(OBJ) 

realclean: clean
	rm -f $(MAIN)

depend:
	(sed '/^# DO NOT DELETE THIS LINE/q' Makefile && \
	${CC} -MM ${SRC} \
	) >Makefile.new
	cp Makefile Makefile.bak
	cp Makefile.new Makefile
	rm -f Makefile.new

# DO NOT DELETE THIS LINE
SegmentBM.o: SegmentBM.c Geometry.h PartitionBoundaryMesh.h \
 Postprocessing.h
Geometry.o: Geometry.c Geometry.h
List.o: List.c List.h
PartitionBoundaryMesh.o: PartitionBoundaryMesh.c Geometry.h List.h
Postprocessing.o: Postprocessing.c Postprocessing.h Geometry.h

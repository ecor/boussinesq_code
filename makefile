################################################################################
# Draft Makefile for THe Boussinesq Equation Model
#################################################################################
# Author: Emanuele Cordano
# Date:   24 jun 2013
# Comment: The binary will be created in the same folder where the makefile is.

# link: http://stackoverflow.com/questions/15965390/simple-makefile-not-taking-include-path
SHELL = /bin/sh
CC=gcc
CFLAGS=-g -Wall
HM = .
BINPATH   = $(HM)/bin
NAME		= boussinesq
EXECUTABLE		= $(BINPATH)/$(NAME)

SOURCEDIR = $(HM)/source


SOURCES= $(SOURCEDIR)/FLUIDTURTLES/alloc.c \
         $(SOURCEDIR)/FLUIDTURTLES/datamanipulation.c \
         $(SOURCEDIR)/FLUIDTURTLES/error.c \
         $(SOURCEDIR)/FLUIDTURTLES/linearalgebra.c \
         $(SOURCEDIR)/FLUIDTURTLES/list.c \
		 $(SOURCEDIR)/FLUIDTURTLES/probability.c \
		 $(SOURCEDIR)/FLUIDTURTLES/random.c \
	     $(SOURCEDIR)/FLUIDTURTLES/statistics.c \
         $(SOURCEDIR)/FLUIDTURTLES/string.c \
		 $(SOURCEDIR)/FLUIDTURTLES/t_io.c \
         $(SOURCEDIR)/FLUIDTURTLES/tensors3D.c \
         $(SOURCEDIR)/FLUIDTURTLES/utilities.c \
         $(SOURCEDIR)/FLUIDTURTLES/write_dem.c \
         $(SOURCEDIR)/ASCII/rw_maps.c \
         $(SOURCEDIR)/ASCII/import_ascii.c \
         $(SOURCEDIR)/ASCII/tabs.c \
         $(SOURCEDIR)/ASCII/write_ascii.c \
         $(SOURCEDIR)/BGEOMETRY/bigcells2.c \
         $(SOURCEDIR)/BGEOMETRY/g_raster2plvector.c \
         $(SOURCEDIR)/BGEOMETRY/geometry.c \
         $(SOURCEDIR)/BGEOMETRY/geometry2.c \
         $(SOURCEDIR)/BGEOMETRY/geometry_attribute.c \
         $(SOURCEDIR)/BGEOMETRY/geometry_freememory.c \
         $(SOURCEDIR)/BGEOMETRY/geometry_io.c \
         $(SOURCEDIR)/BGEOMETRY/geometry_utilities.c \
         $(SOURCEDIR)/BGEOMETRY/sorting.c \
         $(SOURCEDIR)/KeyPalette/additional_read_functions.c \
		 $(SOURCEDIR)/KeyPalette/get_filenames.c \
		 $(SOURCEDIR)/KeyPalette/key.palette.c \
         $(SOURCEDIR)/KeyPalette/read_command_line.c \
         $(SOURCEDIR)/MATH2/boussinesq_matrix.c \
         $(SOURCEDIR)/MATH2/doublevector_utilities.c \
         $(SOURCEDIR)/MATH2/geo_statistic.09375.c \
         $(SOURCEDIR)/MATH2/linear_span.c \
         $(SOURCEDIR)/MATH2/preconditioned_conjugate_gradient.c \
         $(SOURCEDIR)/MATH2/util_math.c \
         $(SOURCEDIR)/Boussinesq/b_readgrid.c \
         $(SOURCEDIR)/Boussinesq/b_solver.c \
         $(SOURCEDIR)/Boussinesq/b_sources.c \
         $(SOURCEDIR)/Boussinesq/b_utilities.c \
         $(SOURCEDIR)/Boussinesq/b_v_advection.c \
         $(SOURCEDIR)/Boussinesq/b_volumes.c \
         $(SOURCEDIR)/Boussinesq/main.c 


INCLUDE= -I$(SOURCEDIR)/FLUIDTURTLES \
         -I$(SOURCEDIR)/ASCII \
         -I$(SOURCEDIR)/BGEOMETRY \
         -I$(SOURCEDIR)/KeyPalette \
         -I$(SOURCEDIR)/MATH2 \
         -I$(SOURCEDIR)/Boussinesq



OBJ=$(SOURCES:.c=.o)

.c.o.: $*.c $*.o 
	$(CC)  $(HEADERS) $($(SOURCEDIR)S) $(CFLAGS) -c $(OBJ) $<  -o $@




.c.o: $*.c $*.h
	$(CC) $(CFLAGS) -c $< $(INCLUDE) -o $@


all: $(NAME)

$(NAME): $(OBJ)
	$(CC) -o $(EXECUTABLE) $(OBJ) -lm 
clean:
	rm -rf *.o *~ $(OBJ)







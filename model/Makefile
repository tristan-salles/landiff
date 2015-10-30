#---------------------------------------------------------------------------
#             AUTHOR TRISTAN SALLES -- UNIVERSITY OF SYDNEY
#---------------------------------------------------------------------------
GFORTRAN    = gfortran
MPI_FORTRAN = mpif90
LD = -lm

# OPTIMISATION
#FFLAGS = -O0 -g  -Wall -fbacktrace -lstdc++ -cpp -fcheck=bounds -finit-real=nan \
	-ffpe-trap=zero,overflow,invalid -ffree-form -fno-common -Wtabs -Wunused-parameter \
	-Wuninitialized  -ffree-line-length-none -fdump-fortran-optimized -fdump-tree-original

FFLAGS = -O3

.SUFFIXES : .o .f90

all: multidiff

multidiff : diffusionParam.o diffusionSolve.o diffusionMain.o
	$(MPI_FORTRAN) $(FFLAGS) -o $@ diffusionMain.o diffusionSolve.o diffusionParam.o

.f90.o :
	$(MPI_FORTRAN) $(FFLAGS) -c $(*F).f90

clean :
	/bin/rm -f *.o multidiff

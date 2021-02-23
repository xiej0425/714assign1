# Type 'make' or 'make game' to compile the serial version.
# Type 'make mpi' to compile the MPI version.
# Type 'make collective' to compile the MPI (with Collective I/O) version.
# Type 'make async' to compile the MPI (with Async I/O) version.
# Type 'make openmp' to compile the MPI + OpenMP version.
# Type 'make cuda' to compile the Cuda version.
CC = gcc
MPICC = mpicc
NVCC = nvcc
FLAGS = -std=c99 -Wall -O3
CPP = g++
MPICXX = mpic++
FLAGSCXX = -std=c++14

1d:
	$(MPICXX) $(FLAGSCXX) ass1_mpi.cpp

2d:
	$(MPICC) $(FLAGS) game_mpi.c -lm

inputgen:
	$(CPP) $(FLAGSCXX) 2d_input_gen.cpp -lm

clean:
	rm -f *.out

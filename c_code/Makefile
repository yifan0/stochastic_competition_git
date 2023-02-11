CC=g++
MPICC=mpic++
CFLAGS=-I. -std=c++11 -fopenmp -O3

BLAS_INC = /usr/local/src/openblas/0.3.12/gcc/Sandy.Bridge/include
BLAS_LIB = /usr/local/src/openblas/0.3.12/gcc/Sandy.Bridge/lib
BCL_BACKEND = -DBCL_BACKEND_MPI

LD_FLAGS = $(GASNET_BACKEND_LIB_FLAGS)

all: sim simple sim1d matrix_basic petsc_sim

sim: sim.cpp
	$(CC) $(CFLAGS) -o sim sim.cpp

simple: simple_sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o simple simple_sim.cpp

sim1d: sim_1d.cpp
	$(CC) $(CFLAGS) -Irandomc -o sim1d sim_1d.cpp

matrix_basic: matrix_basic.cpp
	$(MPICC) -std=gnu++17 -O3 $(BCL_BACKEND) -I../bcl -I$(BLAS_INC) -L$(BLAS_LIB) -lopenblas -o $@ $<

petsc_sim: petsc_sim.c
	$(CC) $(CFLAGS) -Irandomc -L$(PETSC_DIR)/lib -lpetsc -o $@ $< $(LD_FLAGS)

clean:
	rm -f sim simple sim1d matrix_basic petsc_sim


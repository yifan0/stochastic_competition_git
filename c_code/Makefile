CC=g++
MPICC=mpic++
CFLAGS=-I. -std=c++17 -O3

BLAS_INC = /usr/local/src/openblas/0.3.12/gcc/Sandy.Bridge/include
BLAS_LIB = /usr/local/src/openblas/0.3.12/gcc/Sandy.Bridge/lib
BCL_BACKEND = -DBCL_BACKEND_MPI
BOOST_INC := /usr/local/boost/1.71.0/include

LD_FLAGS = $(GASNET_BACKEND_LIB_FLAGS)

all: mpi_sim simple

sim: sim.cpp
	$(CC) $(CFLAGS) -o sim sim.cpp

simple: simple_sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o simple simple_sim.cpp

sim1d: sim_1d.cpp
	$(CC) $(CFLAGS) -Irandomc -o sim1d sim_1d.cpp

bcl_sim: bcl_sim.cpp
	$(MPICC) -std=gnu++17 -O3 $(BCL_BACKEND) -I../bcl -I$(BLAS_INC) -L$(BLAS_LIB) -lopenblas -Irandomc -o $@ $<

mpi_sim: mpi_sim.cpp
	$(MPICC) -std=gnu++17 -O3 -Irandomc -o $@ $<

matrix_basic: matrix_basic.cpp
	$(MPICC) -std=gnu++17 -O3 $(BCL_BACKEND) -I../bcl -I$(BLAS_INC) -L$(BLAS_LIB) -I${BOOST_INC} -lopenblas -o $@ $<

petsc_sim: petsc_sim.c
	$(MPICC) $(CFLAGS) -Irandomc -L$(PETSC_DIR)/lib -lpetsc -o $@ $< $(LD_FLAGS)

dash_sim: dash_narray.cpp
	/home/ekoning2/opt/bin/dash-mpicxx $(CFLAGS) -I/home/ekoning2/scratch/dash/build -L/home/ekoning2/scratch/dash/build/lib -Irandomc -o $@ $<

clean:
	rm -f sim simple sim1d matrix_basic petsc_sim


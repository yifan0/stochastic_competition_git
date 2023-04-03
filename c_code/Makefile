CC=g++
MPICC=mpic++
CFLAGS=-I. -std=c++17 -O3

LD_FLAGS = $(GASNET_BACKEND_LIB_FLAGS)

all: mpi_sim sim

sim: sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o sim sim.cpp

mpi_sim: mpi_sim.cpp
	$(MPICC) -std=gnu++17 -O3 -Irandomc -o $@ $<

clean:
	rm -f sim mpi_sim


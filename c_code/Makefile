CC=g++
MPICC=mpic++
CXXOPTS=-I./include
CFLAGS=-I. -std=c++17 -O3 -Irandomc $(CXXOPTS)

LD_FLAGS = $(GASNET_BACKEND_LIB_FLAGS)
GA_INC=-I/home/ekoning2/scratch/ga/install/include
GA_LINK=-L/home/ekoning2/scratch/ga/build/install/lib -L/home/ekoning2/scratch/ga/build/install/lib64 -lgfortran -lga -lga++ -llapack -lblas -larmci -lblas -lrt -lm

EXEC=mpi_sim tree_sim ga_sim sim

all: $(EXEC)

tree_sim: tree_sim.cpp
	$(CC) $(CFLAGS) -o $@ $<

sim: sim.cpp
	$(CC) $(CFLAGS) -o $@ $<

mpi_sim: mpi_sim.cpp
	$(MPICC) $(CFLAGS) -o $@ $<

ga_sim: ga_sim.cpp
	$(MPICC) $(CFLAGS) -O -g $(GA_INC) -DUSE_MPI -DMPI -o $@ $< $(GA_LINK)

clean:
	rm -f $(EXEC)


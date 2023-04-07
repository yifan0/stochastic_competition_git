CC=g++
MPICC=mpic++
CFLAGS=-I. -std=c++17 -O3

LD_FLAGS = $(GASNET_BACKEND_LIB_FLAGS)
GA_INC=-I/home/ekoning2/scratch/ga/install/include
GA_LINK=-L/home/ekoning2/scratch/ga/build/install/lib -L/home/ekoning2/scratch/ga/build/install/lib64 -lgfortran -lga -lga++ -llapack -lblas -larmci -lblas -lrt -lm

EXEC=mpi_sim sim ga_sim

all: $(EXEC)

sim: sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o sim sim.cpp

mpi_sim: mpi_sim.cpp
	$(MPICC) -std=gnu++17 -O3 -Irandomc -o $@ $<

ga_sim: ga_sim.cpp
	$(MPICC) $(CFLAGS) -O -g $(GA_INC) -Irandomc -DUSE_MPI -DMPI -o $@ $< $(GA_LINK)

clean:
	rm -f $(EXEC)


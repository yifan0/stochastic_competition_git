CC=g++
CFLAGS=-I. -std=c++11 -fopenmp -O3

all: sim simple sim1d

sim: sim.cpp
	$(CC) $(CFLAGS) -o sim sim.cpp

simple: simple_sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o simple simple_sim.cpp

sim1d: sim_1d.cpp
	$(CC) $(CFLAGS) -Irandomc -o sim1d sim_1d.cpp

clean:
	rm sim simple


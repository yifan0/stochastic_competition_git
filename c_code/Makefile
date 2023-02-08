CC=g++
CFLAGS=-I. -std=c++11 -fopenmp -O3

all: sim simple

sim: sim.cpp
	$(CC) $(CFLAGS) -o sim sim.cpp

simple: simple_sim.cpp
	$(CC) $(CFLAGS) -Irandomc -o simple simple_sim.cpp

clean:
	rm sim simple


# Implementations with C++

From `c_code` directory, build by running:

`mkdir build && cd build && cmake .. && make`

## Sequential with Trees

- Code: tree\_sim.cpp
- Executable: tree\_sim

## Sequential with Grid Only

- Code: sim.cpp
- Executable: sim

## Parallel with GlobalArrays

- Code: ga\_sim.cpp
- Executable: ga\_sim

# Sequential Benchmark

- Code: benchmark.cpp
- Executable: benchmark
- Provides timing baseline. Executes steps similar to those of the simulation, but without interleaving them, providing a breakdown of required time.



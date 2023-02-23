#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=run_mpi_%j.log   # Standard output and error log
#SBATCH --partition secondary-eth

module load gcc/7.2.0
module load openmpi/4.1.0-gcc-7.2.0-pmi2

echo "1 nodes with 1 task each"

mpirun ./mpi_sim 512 1 test_mpi_3_


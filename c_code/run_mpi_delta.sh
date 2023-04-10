#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/run_mpi_%j.log   # Standard output and error log
#SBATCH --partition cpu
#SBATCH --account bbkf-delta-cpu

echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./mpi_sim -s 500 -o test_mpi_0_


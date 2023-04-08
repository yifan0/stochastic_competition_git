#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=test_mpi_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction 

SIZE=500

echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./mpi_sim ${SIZE} 1 test_mpi_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}_


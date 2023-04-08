#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=test_ga_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction 

module load openmpi/4.1.0-gcc-7.2.0-pmi2

SIZE=500
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./ga_sim ${SIZE} 1 test_ga_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}_


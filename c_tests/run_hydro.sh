#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --time=0-01:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

SIZE=100
OUTPUT=sim_${SLURM_JOB_ID}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

module load openmpi/4.1.4
cd build && make > null && cd ..

srun time ./build/sim -s ${SIZE} -o ${OUTPUT}


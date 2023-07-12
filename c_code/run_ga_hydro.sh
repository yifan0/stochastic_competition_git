#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=7-00:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

SIZE=1000
OUTPUT=ga_sim_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS_PER_NODE}
echo Tasks = ${SLURM_NTASKS_PER_NODE}
echo Nodes = ${SLURM_JOB_NUM_NODES}

module load openmpi/4.1.4
mpirun ./build/ga_sim -s ${SIZE} -o ${OUTPUT}

module load anaconda3_cpu/4.13.0
python heatmap.py ${OUTPUT}_rep0.csv


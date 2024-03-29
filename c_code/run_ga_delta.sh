#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=121
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition cpu
#SBATCH --account bbkf-delta-cpu

SIZE=500
OUTPUT=test_ga_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./build/ga_sim -s ${SIZE} -o ${OUTPUT}

module load anaconda3_cpu/4.13.0
python heatmap.py ${OUTPUT}_rep0.csv


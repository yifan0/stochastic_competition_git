#!/bin/bash
#SBATCH --job-name=ga_sim1D_hydro    # Job name
#SBATCH --time=2:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

SIZE=262144
OUTPUT=container/1D_262144_4/ga_sim_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS_PER_NODE}
echo Nodes = ${SLURM_JOB_NUM_NODES}

module load openmpi/4.1.4
mpirun ./build/ga_sim1D_hydro -s ${SIZE} -o ${OUTPUT} --specrate=1.0e-4

# module load anaconda3_cpu/4.13.0
# python heatmap.py ${OUTPUT}_rep0.csv


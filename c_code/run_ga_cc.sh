#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction 

module unload anaconda/2022-May/3
module load openmpi/4.1.0-gcc-7.2.0-pmi2

SIZE=500
OUTPUT=test_ga_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./ga_sim -s ${SIZE} -o ${OUTPUT}

module load anaconda/2022-May/3
python reformat_ga_matrix.py ${OUTPUT}_rep0.log
python heatmap.py ${OUTPUT}_rep0.csv


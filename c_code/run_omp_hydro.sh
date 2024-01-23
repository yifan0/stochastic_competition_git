#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=0-07:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

SIZE=100
OUTPUT=ga_sim_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

module load openmpi/4.1.4
module load python
cd build && make > null && cd ..

export OMP_NUM_THREADS=8
srun ./build/omp_sim -s ${SIZE} -o ${OUTPUT} 
python heatmap.py ${OUTPUT}_rep0.csv

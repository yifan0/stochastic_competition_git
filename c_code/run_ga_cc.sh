#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction 

module unload anaconda/2022-May/3
module load openmpi/4.1.0-gcc-7.2.0-pmi2

SIZE=500
OUTPUT=ga_test_for_mixing_${SIZE}_${SLURM_JOB_NUM_NODES}_${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

mpirun ./build/ga_sim -s ${SIZE} -o ${OUTPUT} --specrate=1e-7 

module load anaconda/2022-May/3
for i in {0..99..1}
do
    python heatmap.py ${OUTPUT}_rep0_checkpoint${i}.csv
done

python heatmap.py ${OUTPUT}_rep0.csv


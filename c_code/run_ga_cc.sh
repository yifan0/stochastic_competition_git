#!/bin/bash
#SBATCH --job-name=ga_sim    # Job name
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --output=logs/test_ga_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction 
#SBATCH --account wgropp-ic

module load openmpi/5.0.1-gcc-13.3.0

SIZE=500
DIMS=1
OUTPUT=ga_test_for_mixing_size${SIZE}_${DIMS}D_nodes${SLURM_JOB_NUM_NODES}_tasks${SLURM_NTASKS}
echo Tasks = ${SLURM_NTASKS}
echo Nodes = ${SLURM_JOB_NUM_NODES}

cd /u/ekoning2/stochastic_competition_git/c_code/
# make
mpirun -n ${SLURM_NTASKS} ./build_test/ga_sim -s ${SIZE} -o ${OUTPUT} --specrate=1e-7 -d ${DIMS} -f -u

#module load anaconda/2022-May/3
#for i in {0..99..1}
#do
#    python heatmap.py ${OUTPUT}_rep0_checkpoint${i}.csv
#done

#python heatmap.py ${OUTPUT}_rep0.csv


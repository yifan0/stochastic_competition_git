#!/bin/bash
#SBATCH --job-name=omp_sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/omp_sim_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

export OMP_NUM_THREADS=10
SIZE=100
OUTPUT=test_omp_${SIZE}
module unload anaconda/2022-May/3
cd build && make && cd ..
/usr/bin/time -v ./build/omp_sim -s ${SIZE} -o ${OUTPUT}
module load anaconda/2022-May/3
python heatmap.py ${OUTPUT}_rep0.csv


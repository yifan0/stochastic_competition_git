#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/sim_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

SIZE=500
OUTPUT=test_sequential_${SIZE}
module unload anaconda/2022-May/3
cd build && make && cd ..
/usr/bin/time -v ./build/sim -s ${SIZE} -o ${OUTPUT}
module load anaconda/2022-May/3
python heatmap.py ${OUTPUT}_rep0.csv


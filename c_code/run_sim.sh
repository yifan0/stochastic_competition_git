#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/sim_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

SIZE=900

/usr/bin/time -v ./build/sim -s ${SIZE} -o test_sequential_${SIZE}_


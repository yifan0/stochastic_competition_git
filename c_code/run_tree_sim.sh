#!/bin/bash
#SBATCH --job-name=tree    # Job name
#SBATCH --ntasks=1
#SBATCH --time=01:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/tree_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

SIZE=100
./build/tree_sim -s ${SIZE} -o test_sequential_tree_${SIZE}


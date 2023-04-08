#!/bin/bash
#SBATCH --job-name=simple    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=simple_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

SIZE=1000

/usr/bin/time -v ./sim ${SIZE} 1 test_sequential_${SIZE}


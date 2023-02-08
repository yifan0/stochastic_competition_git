#!/bin/bash
#SBATCH --job-name=run_sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=12:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=run_sim_%j.log   # Standard output and error log
#SBATCH --partition cs

make clean && make
/usr/bin/time -v /home/ekoning2/stochastic_competition_git/c_code/sim 512 1 test10.txt


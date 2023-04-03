#!/bin/bash
#SBATCH --job-name=simple    # Job name
#SBATCH --ntasks=1
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=simple_%j.log   # Standard output and error log
#SBATCH --partition secondary

/usr/bin/time -v /home/ekoning2/stochastic_competition_git/c_code/simple 1024 1 testI.txt


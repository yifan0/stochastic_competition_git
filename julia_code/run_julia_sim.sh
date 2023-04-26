#!/bin/bash
#SBATCH --job-name=run_sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=00:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --output=run_sim_%j.log   # Standard output and error log
#SBATCH --partition cs

/usr/bin/time -v /home/ekoning2/julia/julia-1.8.5/bin/julia sim_git.jl


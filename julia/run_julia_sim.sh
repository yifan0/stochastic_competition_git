#!/bin/bash
#SBATCH --job-name=run_sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=00:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --output=sim_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

module load julia/1.8.5

/usr/bin/time -v julia -t 16 sim.jl


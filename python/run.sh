#!/bin/bash
#SBATCH --job-name=mesa_sim    # Job name
#SBATCH --time=0-05:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=sim_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

time python sim.py


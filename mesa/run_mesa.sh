#!/bin/bash
#SBATCH --job-name=mesa_sim    # Job name
#SBATCH --time=0-05:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=test_mesa_%j.log   # Standard output and error log
#SBATCH --partition sandybridge

time python run.py


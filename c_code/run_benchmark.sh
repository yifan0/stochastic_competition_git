#!/bin/bash
#SBATCH --job-name=sim    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/benchmark_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

./build/benchmark


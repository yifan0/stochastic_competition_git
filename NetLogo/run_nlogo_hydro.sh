#!/bin/bash
#SBATCH --job-name=nlogo_sim    # Job name
#SBATCH --time=0-05:20:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition sandybridge

time ~/NetLogo/NetLogo-6.3.0/NetLogo_Console --headless --model sim.nlogo --experiment experiment --table table.csv


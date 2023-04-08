#!/bin/bash
#SBATCH --job-name=wrap    # Job name
#SBATCH --ntasks=1
#SBATCH --time=04:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=wrap_%j.log   # Standard output and error log
#SBATCH --partition eng-instruction

for SIZE in {100..1000..100}
do
    /home/ekoning2/stochastic_competition_git/c_code/wrap_sim ${SIZE} 1 test_sequential_${SIZE}
done


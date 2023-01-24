#!/bin/bash
#SBATCH --job-name=simple    # Job name
#SBATCH --ntasks=1
#SBATCH --time=02:00:00                 # Time limit hrs:min:sec
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=simple_%j.log   # Standard output and error log
#SBATCH --partition cs

module load boost/1.71.0
module load intel/18.0

cat /home/ekoning2/stochastic_competition_git/c_code/simple_sim.cpp

export OMP_NUM_THREADS=1
echo "Threads = " $OMP_NUM_THREADS
make clean && make
/usr/bin/time -v /home/ekoning2/stochastic_competition_git/c_code/simple 512 1 test9.txt


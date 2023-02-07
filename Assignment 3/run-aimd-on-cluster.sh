#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                    # Project
#SBATCH -J AIMD                            # Name of the job
#SBATCH -N 1                               # Use 1 node
#SBATCH -n 32                              # Use all cores on that node
#SBATCH -t 00:01:00                        # Maximum time
#SBATCH -o "./Assignment 3/logs/std.out"   # stdout goes to this file
#SBATCH -e "./Assignment 3/logs/err.out"   # stderr goes to this file

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

mpiexec gpaw python 'Assignment 3/run-aimd.py'

# For 3 simulation timesteps
# 1 core:   3:39
# 8 cores:  0:39
# 32 cores: 0:18
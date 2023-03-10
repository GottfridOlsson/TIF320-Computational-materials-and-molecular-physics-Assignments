#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                    # Project
#SBATCH -J A4_T6                           # Name of the job
#SBATCH -N 1                               # Use 1 node
#SBATCH -n 32                              # Use 32 cores on that node
#SBATCH -t 10:00:00                        # Maximum time
#SBATCH -o "Assignment 4/logs/T6_std.out"  # stdout goes to this file
#SBATCH -e "Assignment 4/logs/T6_err.out"  # stderr goes to this file

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

mpiexec gpaw python 'Assignment 4/task6.py'
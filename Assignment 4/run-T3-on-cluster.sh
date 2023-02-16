#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                    # Project
#SBATCH -J A4_T3                           # Name of the job
#SBATCH -N 1                               # Use 1 node
#SBATCH -n 2                               # Use 2 cores on that node
#SBATCH -t 02:00:00                        # Maximum time
#SBATCH -o "Assignment 4/logs/T3_std.out"   # stdout goes to this file
#SBATCH -e "Assignment 4/logs/T3_err.out"   # stderr goes to this file

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

python3 'Assignment 4/task3.py'


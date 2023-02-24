#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                    # Project
#SBATCH -J A4_T5                           # Name of the job
#SBATCH -N 1                               # Use 1 node
#SBATCH -n 4                               # Use 16 cores on that node
#SBATCH -t 00:30:00                        # Maximum time
#SBATCH -o "logs/T5_std.out"  # stdout goes to this file
#SBATCH -e "logs/T5_err.out"  # stderr goes to this file

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

mpiexec gpaw python 'task5.py'
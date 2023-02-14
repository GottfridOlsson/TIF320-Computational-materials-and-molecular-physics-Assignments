#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                    # Project
#SBATCH -J T1_lattice_params               # Name of the job
#SBATCH -N 1                               # Use 1 node
#SBATCH -n 1                               # Use 1 core on that node
#SBATCH -t 01:00:00                        # Maximum time
#SBATCH -o "./Assignment 4/logs/std.out"   # stdout goes to this file
#SBATCH -e "./Assignment 4/logs/err.out"   # stderr goes to this file

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

python3 'Assignment 4/task1.py'
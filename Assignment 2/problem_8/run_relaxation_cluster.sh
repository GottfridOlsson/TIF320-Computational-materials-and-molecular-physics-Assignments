#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2                         # Project
#SBATCH -J GPAW_P8                              # Name of the job
#SBATCH -N 1                                    # Use 1 node
#SBATCH -n 1                                    # Use only 1 core on that node
#SBATCH -t 01:00:00                             # Maximum time
#SBATCH -o "./Assignment 2/problem_8/std.out"   # stdout goes to this file
#SBATCH -e "./Assignment 2/problem_8/err.out"   # stderr goes to this file
#SBATCH --array 0-11                            # Which settings to use

# Run from repository root using 
# sbatch "./Assignment 2/problem_8/run_relaxation_cluster.sh"

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

setting_no=$SLURM_ARRAY_TASK_ID

# Run relaxation
python3 "Assignment 2/problem_8/relaxation.py" $setting_no
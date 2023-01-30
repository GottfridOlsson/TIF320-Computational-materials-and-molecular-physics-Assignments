#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2     # Project
#SBATCH -J TIF320-A2-P1     # Name of the job
#SBATCH -N 1                # Use 1 node
#SBATCH -n 1                # Use only 1 core on that node
#SBATCH -t 10:00:00         # Maximum time
#SBATCH -o std.out          # stdout goes to this file
#SBATCH -e err.out          # stderr goes to this file
#SBATCH --array 6-8         # Which number of atoms to simulate for

# Run using 
# sbatch run_ga_cluster.sh

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

natoms=$SLURM_ARRAY_TASK_ID

echo "Running genetic algorithm for natoms=$natoms."

# Go to folder corresponding to this atom number
folder="natoms_$natoms"
cd $folder

# Run genetic algorithm
python3 ../../ga.py

echo "Simulation done for natoms=$natoms!"

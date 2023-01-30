#!/usr/bin/env bash

# Load modules
module purge
module load SciPy-bundle/2022.05-foss-2022a
module load ASE/3.22.1-foss-2022a

# Run program with 6, 7 and 8 atoms in different folders
for natoms in 6 7 8
do
    echo "Running initialization for $natoms atoms..."
    folder="natoms_$natoms"
    mkdir $folder
    cd $folder
    python3 ../../initialization.py $natoms
    cd ..
done
echo "Done!"
#!usr/bin/env bash

for natoms in 6 7 8
do
    db_path=../problem_1/natoms_${natoms}/gadb.db
    out_path=./databases/natoms_${natoms}.txt
    ase db $db_path -c++ -L 0 -s=energy > $out_path
done
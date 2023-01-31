from ase.db import connect
from ase.io import write

for natoms in [6, 7, 8]:

    assignment_path = './TIF320-Computational-materials-and-molecular-physics-Assignments/Assignment 2'
    db_path = assignment_path + f'/problem_1/natoms_{natoms}/gadb.db'
    save_path = assignment_path + f'/problem_5/structures/natoms_{natoms}'

    # Connect to database
    db = connect(db_path)

    # Create a iterator that returns configurations sorted by energy
    atoms_iterator_gpaw = db.select("calculator=gpaw", sort="energy")
    atoms_iterator_lj = db.select("calculator=lennardjones", sort="energy")

    # Get the first value of the iterator as groundstate
    atoms_groundstate_gpaw =  next(atoms_iterator_gpaw).toatoms()
    atoms_groundstate_lj =  next(atoms_iterator_lj).toatoms()

    # Save the result
    write(save_path + '_groundstate_gpaw.xyz', atoms_groundstate_gpaw)
    write(save_path + '_groundstate_gpaw.png', atoms_groundstate_gpaw)
    write(save_path + '_groundstate_lennardjones.xyz', atoms_groundstate_lj)
    write(save_path + '_groundstate_lennardjones.png', atoms_groundstate_lj)
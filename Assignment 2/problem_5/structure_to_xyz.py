from ase.db import connect
from ase.io import write

for natoms in [6, 7, 8]:

    # Connect to database
    db = connect(f'./TIF320-Computational-materials-and-molecular-physics-Assignments/Assignment 2/problem_1/natoms_{natoms}/gadb.db')

    # Create a iterator that returns configurations sorted by energy
    atoms_iterator = db.select(sort="energy")

    # Get the first value of the iterator
    atoms =  next(atoms_iterator).toatoms()

    # Save the result
    write(f'./TIF320-Computational-materials-and-molecular-physics-Assignments/Assignment 2/problem_5/structures/natoms_{natoms}.xyz', atoms)
    write(f'./TIF320-Computational-materials-and-molecular-physics-Assignments/Assignment 2/problem_5/structures/natoms_{natoms}.png', atoms)
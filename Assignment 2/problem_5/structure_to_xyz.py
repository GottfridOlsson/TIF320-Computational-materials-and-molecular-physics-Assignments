from ase.db import connect
from ase.io import write

db = connect('../problem_1/natoms_6/gadb.db')
atoms = db.get('id=79').toatoms ()
write('./structures/natoms_6.xyz', atoms)
print(atoms.structure)
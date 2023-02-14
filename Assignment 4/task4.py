# Construct CO molecule:
# https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/wavefunctions/wavefunctions/wavefunctions.html
from ase.build import molecule
CO = molecule('CO')

cell_length = 12 #angstrom

atom = Atoms('O', cell=[6, 6, 6], pbc=False)
atom.center()

grid_spacing = 0.2 #anstrom


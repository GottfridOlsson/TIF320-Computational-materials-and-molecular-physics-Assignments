# Construct surface
# Example from https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/surface/surface.html
from ase.visualize import view
from ase.build import fcc100
s = fcc100('Al', (1, 1, 5))
view(s, repeat=(4, 4, 1))
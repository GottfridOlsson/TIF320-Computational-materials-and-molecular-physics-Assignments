# Adapted from 
# https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/aluminium/aluminium.html
# Bulk Al(fcc) test
from ase import Atoms
from ase.io import write
from gpaw import GPAW, PW

name = 'Al-fcc'
a = 4.05  # initial fcc lattice parameter
b = a / 2

bulk = Atoms('Al',
             cell=[[0, b, b],
                   [b, 0, b],
                   [b, b, 0]],
             pbc=True)

write("bulk_Al.png", bulk)

k = 4
calc = GPAW(mode=PW(300),       # cutoff
            kpts=(k, k, k),     # k-points
            txt=name + '.txt')  # output file

bulk.calc = calc

energy = bulk.get_potential_energy()
calc.write(name + '.gpw')
print('Energy:', energy, 'eV')
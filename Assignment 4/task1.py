# Adapted from 
# https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/aluminium/aluminium.html

from ase import Atoms
from ase.io import write
from gpaw import GPAW, PW
import util


# Elements and FCC lattice parameters to explore
elements = ['Au', 'Pt', 'Rh']
a_guesses = [3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4]

# Number of k-points and cutoff energy for plane waves
k = 12
cutoff_energy = 450

for element in elements:
      #print(f"Calculating for element {element}")
      energies = []
      for a in a_guesses:
            #print(f"Lattice parameter a={a}")

            file_name = f"Assignment 4/output/{element}_calculation_a={a:.1f}"

            # Create bulk FCC with periodic boundary conditions (pbc)
            b = a/2
            bulk = Atoms(element,
                        cell=[[0, b, b],
                              [b, 0, b],
                              [b, b, 0]],
                        pbc=True)

            # Set calculator
            calc = GPAW (xc   = 'PBE',               # exchange-correlation
                        mode = PW(cutoff_energy),    # cutoff
                        kpts = (k, k, k),            # k-points
                        txt  = f'{file_name}.txt')    # output file

            # Calculate potential energy of structure
            bulk.calc = calc
            energy = bulk.get_potential_energy()
            calc.write(f'{file_name}.gpw')
            energies.append(energy)
      
      util.print_arrays_to_CSV(f"Assignment 4/output/TIF320_T1_{element}_energy_vs_lattice_parameter.csv", 
                                "Lattice parameter a (Å)", a_guesses, 
                                "Potential energy for {element} (eV)", energies, 
                                print_message=True)
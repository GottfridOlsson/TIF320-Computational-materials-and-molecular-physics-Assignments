# Adapted from 
# https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/aluminium/aluminium.html

from ase import Atoms
from ase.io import write
from gpaw import GPAW, PW
import util
import numpy as np


# Elements and FCC lattice parameters to explore
elements = ['Au', 'Pt', 'Rh']
a_theoretical = [4.0782, 3.9236, 3.8032] # angstrom
a_step = 0.01

# Number of k-points and cutoff energy for plane waves (from problem description)
k = 12
cutoff_energy = 450

for i, element in enumerate(elements):
      #print(f"Calculating for element {element}")
      energies = []
      a_trials = np.arange(a_theoretical[i]-4*a_step, a_theoretical[i]+5*a_step, step=a_step)

      for a in a_trials:
            #print(f"Lattice parameter a={a}")

            file_name = f"Assignment 4/output/{element}_calculation_a={a:.2f}"

            # Create bulk FCC with periodic boundary conditions (pbc)
            b = a/2
            bulk = Atoms(element,
                        cell=[[0, b, b],
                              [b, 0, b],
                              [b, b, 0]],
                        pbc=True)

            #OR MORE SIMPLY: "bulk = bulk(element, 'fcc', a=a)" with "from ase.build import bulk"

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
      
      util.print_arrays_to_CSV(f"Assignment 4/output/TIF320_A4_T1_{element}_energy_vs_lattice_parameter_step0.01.csv", 
                                f"Lattice parameter for {element} (Å)", a_trials, 
                                f"Potential energy for {element} (eV)", energies, 
                                print_message=True)
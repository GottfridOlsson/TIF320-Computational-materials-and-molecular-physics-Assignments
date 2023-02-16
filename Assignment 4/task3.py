# Construct surface
# Example from https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/surface/surface.html
# also: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

from ase.build import fcc111
from ase.io import write
from gpaw import GPAW, PW
from ase.optimize import GPMin
import util


elements = ['Au', 'Pt', 'Rh']
a_from_T1 = [4.177, 3.970, 3.840] #angstrom
output_path_start = 'Assignment 4/output_T3/'
potential_energy_of_surface = []

for i, element in enumerate(elements):
    print(f"i={i}, element={element}")

    # Construct (111)-surface for Au, Pt, Rh
    # 3-layered 3X3 surface cell, 6 angstrom of vacuum in +-z-direction
    surface = fcc111(element, (3, 3, 3), a=a_from_T1[i], vacuum=6.0)
    write(f"{output_path_start}{element}_111-surface_initialized.png", surface)

    # Calculator
    k = (4, 4, 1) #k-sampling according to problem description
    cutoff_energy = 450 # eV
    calc = GPAW(mode=PW(cutoff_energy),
                kpts=k,
                txt=f"{output_path_start}{element}_111-surface_calc.txt")

    # Relax surface (code and settings from Assignment 2, 'relaxation.py')
    surface.set_calculator(calc)
    dyn = GPMin(surface, 
                trajectory=f"{output_path_start}GPMin_{element}_111-surface.traj", 
                logfile=f"{output_path_start}GPMin_{element}_111-surface.log")
    dyn.run(fmax=0.001, steps=100) #0.001, 100
    write(f"{output_path_start}{element}_111-surface_relaxed.png", surface)
    calc.write(f"{output_path_start}{element}_111-surface_relaxed.gpaw")
    
    # Get surface energy (=potential energy?)
    potential_energy_of_surface.append(surface.get_potential_energy())


util.print_arrays_to_CSV("Assignment 4/output_T3/TIF320_A4_T3_surface_energy_Au-Pt-Rh.csv", 
                            "Element symbol", elements, 
                            "Potential energy of surface (eV)", potential_energy_of_surface, 
                            print_message=True)
    


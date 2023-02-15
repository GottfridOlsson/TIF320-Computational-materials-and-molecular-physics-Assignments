# Construct surface
# Example from https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/surface/surface.html

from ase.build import fcc111
from ase.io import write
from gpaw import GPAW, PW

elements = ['Au', 'Pt', 'Rh']
a_from_T1 = [4.177, 3.970, 3.840] #angstrom
output_path_start = 'Assignment 4/output_T2/'



for i, element in enumerate(elements):

    # Construct (111)-surface for Au, Pt, Rh
    # 3-layered 3X3 surface cell, 6 angstrom of vacuum in +-z-direction
    surface = fcc111(element, (3, 3, 3), a=a_from_T1[i], vacuum=6.0)
    write(f"{output_path_start}{element}_111-surface_test.png", surface)

    # Calculator
    k = (4, 4, 1) #k-sampling according to problem description
    cutoff_energy = 450 # eV
    calc = GPAW(mode=PW(cutoff_energy), kpts=k, txt=f"{output_path_start}{element}_111-surface_calc.txt")

    # Relax surfaces, determine the surface energy for each metal
    # TODO
    


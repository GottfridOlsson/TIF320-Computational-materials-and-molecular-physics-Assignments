# Construct CO molecule:

"""
Calculate the energy of the O2 molecule and the CO molecule.
Your k-point sampling should consist only of the gamma-point.

Hint 2: Are the molecules spin-polarized?
"""
from ase.build import molecule
from ase.io import write
from gpaw import GPAW, PW
import util

molecule_names = ['CO', 'O2']
cell_length = 12 #angstrom
cutoff_energy = 450 #electron volt
output_path_start = "Assignment 4/output_T4/"
potential_energy = []

for molecule_name in molecule_names:

    # Create molecule as atoms object, center, and set calculator
    atoms = molecule(molecule_name, 
                    cell=(cell_length, cell_length, cell_length))
    atoms.center()
    write(f"{output_path_start}{molecule_name}_initial_molecule_structure.png", atoms, rotation='10z,-80x')
    
    calc = GPAW(mode=PW(cutoff_energy),
                kpts={'gamma': True},
                txt=f"{output_path_start}{molecule_name}_calc.txt")

    atoms.set_calculator(calc)
    potential_energy.append(atoms.get_potential_energy())

    # Spin polarized: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html (GPAW will take care of that?)

util.print_arrays_to_CSV(f"{output_path_start}TIF320_A4_T4_potential-energy_of_CO_O2_gas.csv", 
                        "Molecule symbol", molecule_names, 
                        "Potential energy (eV)", potential_energy, 
                        print_message=True)



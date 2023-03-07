from ase.build import molecule
from ase.io import write
from gpaw import GPAW, PW
from ase.optimize import GPMin
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
    
    calc = GPAW(xc='PBE',
                mode=PW(cutoff_energy),
                kpts={'gamma': True},
                txt=f"{output_path_start}{molecule_name}_calc.txt")

    atoms.set_calculator(calc)
    dyn = GPMin(atoms, 
                trajectory=f"{output_path_start}GPMin_{molecule_name}.traj", 
                logfile=f"{output_path_start}GPMin_{molecule_name}.log")
    dyn.run(fmax=0.01, steps=25) #0.001, 100

    potential_energy.append(atoms.get_potential_energy())
    write(f"{output_path_start}{molecule_name}_relaxed_molecule_structure.xyz", atoms)

    util.print_arrays_to_CSV(f"{output_path_start}TIF320_A4_T4_potential-energy_of_CO_O2_gas_relaxed.csv", 
                            "Molecule symbol", molecule_names, 
                            "Potential energy (eV)", potential_energy, 
                            print_message=True)
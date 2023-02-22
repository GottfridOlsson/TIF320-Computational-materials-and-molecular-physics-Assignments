# Inspired by:
# [1]: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

from gpaw import GPAW, PW
from ase import Atoms
from ase.io import read, write
from ase.build import add_adsorbate, molecule, fcc111
from ase.optimize import GPMin
from ase.constraints import FixAtoms


surface   = "111-surface"
surface_names = ["Au", "Pt", "Rh"]
adsorbate_names = ["CO", "O"]

a_T1 = [4.177, 3.970, 3.840] # angstrom, from Task 1 (2023-02-22)
z_distance_adsorbant = 2 # angstrom, from the problem description
positions = ['ontop', 'bridge', 'fcc', 'hcp']
output_path_start = "Assignment 4/output_T6/"

for i, surface_name in enumerate(surface_names):

    a = a_T1[i]

    for adsorbate_name in adsorbate_names:
        
        # Build adsorbate
        if adsorbate_name == "CO": adsorbate = molecule(adsorbate_name)
        if adsorbate_name == "O":  adsorbate = Atoms(adsorbate_name)

        for position in positions:
            
            #print(f"{surface_name}, {a} Å, {adsorbate_name}, {position}")
            # Build surfaces and add adsorbate
            surface = fcc111(surface_name, a=a, size=[3,3,3], vacuum=6.0)
            add_adsorbate(surface, adsorbate, height=z_distance_adsorbant, position=position, mol_index=-1)

            # Constrain all atoms except the adsorbate: (from source [1])
            fixed = list(range(len(surface) - 1))
            surface.constraints = [FixAtoms(indices=fixed)]

            
            # Set calculator 
            calculator = GPAW(xc='PBE',
                              mode=PW(450),
                              kpts=(4, 4, 1), #kpts taken from Task 3, this is a guess
                              txt=f"{output_path_start}GPAW_{surface_name}_{adsorbate_name}_{position}.txt")
            surface.set_calculator(calculator)

            # Relax until the forces acting on each atom are all below 0.1 eV/Å
            dyn = GPMin(surface, 
                        trajectory=f"{output_path_start}GPMin_{surface_name}_{adsorbate_name}_{position}.traj", 
                        logfile=f"{output_path_start}GPMin_{surface_name}_{adsorbate_name}_{position}.log")
            dyn.run(fmax=0.1, steps=100)
            #E_pot = i
            E_pot = surface.get_potential_energy() 
            
            # Output
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed.xyz", surface)
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed.png", surface)
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed_angled-view.png", surface,  rotation='10z,-80x, 5y')
            
            with open(f"{output_path_start}E_pot_surface_and_adsorbant.txt", 'a') as file:
                file.write(f"{surface_name}, {adsorbate_name}, {position}, {E_pot}\n")


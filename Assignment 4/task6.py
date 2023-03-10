# Inspired by:
# [1]: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html


from gpaw import GPAW, PW
from ase import Atoms
from ase.io import read, write
from ase.build import add_adsorbate, molecule, fcc111
from ase.optimize import GPMin
from ase.constraints import FixAtoms
from ase.parallel import paropen


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

            do_calculation = 0

            #if surface_name=='Au' and adsorbate_name=='CO' and position == 'ontop': do_calculation = 1
            #if surface_name=='Pt' and adsorbate_name=='CO' and position == 'ontop': do_calculation = 1
            #if surface_name=='Rh' and adsorbate_name=='CO' and position == 'ontop': do_calculation = 1

            #if surface_name=='Au' and adsorbate_name=='O' and position == 'fcc': do_calculation = 1
            if surface_name=='Pt' and adsorbate_name=='O' and position == 'fcc': do_calculation = 1
            #if surface_name=='Rh' and adsorbate_name=='O' and position == 'fcc': do_calculation = 1
            
            if not do_calculation:
                continue
          
            print(f"{surface_name}, {a} Å, {adsorbate_name}, {position}")
            # Build surfaces and add adsorbate
            
            surface = fcc111(surface_name, a=a, size=[3,3,3], vacuum=6.0)
            add_adsorbate(surface, adsorbate, height=z_distance_adsorbant, position=position, mol_index=-1)

            # Constrain all atoms except the adsorbate: (from source [1])
            #fixed = list(range(len(surface) - 1))
            #surface.constraints = [FixAtoms(indices=fixed)]

            
            # Set calculator 
            calculator = GPAW(xc='PBE',
                              mode=PW(450),
                              kpts=(4, 4, 1), #kpts taken from Task 3, this is a guess
                              txt=f"{output_path_start}GPAW_{surface_name}_{adsorbate_name}_{position}_notFixedAtoms.txt")
            surface.set_calculator(calculator)

            # Relax until the forces acting on each atom are all below 0.1 eV/Å
            dyn = GPMin(surface, 
                        trajectory=f"{output_path_start}GPMin_{surface_name}_{adsorbate_name}_{position}_notFixedAtoms.traj", 
                        logfile=f"{output_path_start}GPMin_{surface_name}_{adsorbate_name}_{position}_notFixedAtoms.log")
            
            # UPDATE after feedback: converge O on Au better for fcc
            # let fmax=0.01 for all cases
            dyn.run(fmax=0.01, steps=30)
            
            E_pot = surface.get_potential_energy() 
            
            # Output
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed_notFixedAtoms.xyz", surface)
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed_notFixedAtoms.png", surface)
            write(f"{output_path_start}{surface_name}_{adsorbate_name}_{position}_relaxed_notFixedAtoms_angled-view_.png", surface,  rotation='10z,-80x, 5y')
            
            print(f"{surface_name}, {a} Å, {adsorbate_name}, {position}, E_pot={E_pot}")
            # Paropen only writes to the file for when world.rank==0, i.e. for the process that gets highest rank (only writes once to file, and not once per core used in calculation)
            with paropen(f"{output_path_start}E_pot_surface_and_adsorbant.txt", 'a') as file:
                file.write(f"{surface_name}, {adsorbate_name}, {position}, {E_pot}\n")
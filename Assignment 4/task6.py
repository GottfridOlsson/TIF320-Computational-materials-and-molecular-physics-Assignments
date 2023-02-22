# Inspired by:
# [1]: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

from gpaw import GPAW, PW
from ase.io import read, write
from ase.build import add_adsorbate
from ase.optimize import GPMin
from ase.constraints import FixAtoms


# Read relaxed surfaces from Task 3 and relaxed molecules from Task 4
Au_111surface = read('Assignment 4/output_T3/Au_111-surface_relaxed.xyz')
Pt_111surface = read('Assignment 4/output_T3/Pt_111-surface_relaxed.xyz')
Rh_111surface = read('Assignment 4/output_T3/Rh_111-surface_relaxed.xyz')
CO_moleule    = read('Assignment 4/output_T4/CO_relaxed_molecule_structure.xyz')
O_atom        = atom() #TODO: create an O atom

surfaces   = [Au_111surface, Pt_111surface, Rh_111surface]
adsorbates = [CO_moleule, O_atom]

z_distance_adsorbant = 2 # angstrom, from the problem description
positions = ['top', 'hcp', 'bridge', 'fcc']
output_path_start = "Assignment 4/output_T6/"

for surface in surfaces:
    for adsorbate in adsorbates:
        for position in positions:

            add_adsorbate(surface, adsorbate, height=z_distance_adsorbant, position=position)

            # Constrain all atoms except the adsorbate: (from source [1])
            fixed = list(range(len(surface) - 1))
            surface.constraints = [FixAtoms(indices=fixed)]

            # Set calculator 
            calculator = GPAW(xc='PBE',
                              mode=PW(450),
                              kpts=(4, 4, 1), #kpts taken from Task 3, this is a guess
                              txt=f"{output_path_start}GPAW_{surface}_{adsorbate}_{position}.txt")
            surface.set_calculator(calculator)

            # Relax until the forces acting on each atom are all below 0.1 eV/Å
            dyn = GPMin(surface, 
                        trajectory=f"{output_path_start}GPMin_{surface}_{adsorbate}_{position}.traj", 
                        logfile=f"{output_path_start}GPMin_{surface}_{adsorbate}_{position}.log")
            dyn.run(fmax=0.1, steps=100)

            E_pot = surface.get_potential_energy()
            
            # save energy of adbsorbed CO and O
            

            # Output
            write(f"{output_path_start}{surface}_{adsorbate}_{position}_relaxed.xyz", surface)
            write(f"{output_path_start}{surface}_{adsorbate}_{position}_relaxed.png", surface)
            write(f"{output_path_start}{surface}_{adsorbate}_{position}_relaxed_angled-view.png", surface,  rotation='10z,-80x, 5y')
            # TODO:
            # calculate and print the E_ads (see problem description) to .txt


            #atoms.calc = EMT()
            #opt = BFGS(atoms, logfile=None)
            #opt.run(fmax=0.01)

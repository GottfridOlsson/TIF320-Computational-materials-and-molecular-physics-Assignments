from ase.io import read, write
from gpaw import GPAW, FermiDirac
from ase.optimize import GPMin
import sys

def run_relaxation(xyz_path, save_path, calc):
    
    # Read groundstate and excited state
    atoms = read(xyz_path)

    # Optimize
    atoms.set_calculator(calc)
    dyn = GPMin(
        atoms, 
        trajectory=save_path + 'relax_ref.traj', 
        logfile=save_path + 'relax_ref.log')
    dyn.run(fmax=0.02, steps=100)

    write(save_path + 'relaxed_configuration.xyz', atoms)
    write(save_path + 'relaxed_configuration.png', atoms)
    calc.write(save_path + 'calc_state.gpw')



if __name__ == "__main__":

    # Parse input
    argc = len(sys.argv)
    if (argc != 1): print("Usage: python3 run_task8.py setting_number")
    setting_number = int(sys.argv[0])


    groundstate_path  = './Assignment 2/Na6-structures/christmas-tree.xyz'
    excited_path = './Assignment 2/Na6-structures/half-decahedron.xyz'

    if setting_number == 0:
        save_path = './Assignment 2/problem_8/LCAO/groundstate/'
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands=10,
                h=0.25,
                txt=save_path + 'gpaw-out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp'
            )
        )

    elif setting_number == 1:
        save_path = './Assignment 2/problem_8/LCAO/groundstate/'
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands=10,
                h=0.25,
                txt=save_path + 'gpaw-out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp'
            )
        )

    elif setting_number == 1:
        save_path = './Assignment 2/problem_8/LCAO/groundstate/'
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands=10,
                h=0.25,
                txt=save_path + 'gpaw-out.txt',
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode='lcao',
                basis='dzp'
            )
        )
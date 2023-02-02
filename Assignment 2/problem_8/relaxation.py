from ase.io import read, write
from gpaw import GPAW, FermiDirac, PW
from ase.optimize import GPMin
import sys
import os
import time

def run_relaxation(xyz_path, save_path, calc):
    
    # Read groundstate and excited state
    atoms = read(xyz_path)

    # Optimize
    atoms.set_calculator(calc)
    dyn = GPMin(
        atoms, 
        trajectory=save_path + 'relax_ref.traj', 
        logfile=save_path + 'relax_ref.log')
    dyn.run(fmax=0.001, steps=100)

    write(save_path + 'relaxed_configuration.xyz', atoms)
    write(save_path + 'relaxed_configuration.png', atoms)
    calc.write(save_path + 'calc_state.gpw')


if __name__ == "__main__":

    # Parse input
    argc = len(sys.argv)
    if (argc != 2): 
        print("Usage: python3 run_task8.py setting_number")
        exit()
    setting_number = int(sys.argv[1])

    # Define paths to where to find the initial structures
    groundstate_path  = 'Assignment 2/Na6-structures/christmas-tree.xyz'
    excited_path = 'Assignment 2/Na6-structures/half-decahedron.xyz'

    if setting_number == 0:
        name = "LCAO/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'dzp'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    
    if setting_number == 1:
        name = "LCAO/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'dzp'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    
    if setting_number == 2:
        name = "planewave_50/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = PW(50)
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    
    if setting_number == 3:
        name = "planewave_50/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = PW(50)
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")


    if setting_number == 4:
        name = "planewave_500/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = PW(500)
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    
    if setting_number == 5:
        name = "planewave_500/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = PW(500)
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")

    if setting_number == 6:
        name = "finitedifference/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'fd'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    

    if setting_number == 7:
        name = "finitedifference/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'fd'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")

    if setting_number == 8:
        name = "LCAO_basis_szp/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 6,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'sz(dzp)'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    

    if setting_number == 9:
        name = "LCAO_basis_szp/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 6,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'sz(dzp)'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")

    if setting_number == 10:
        name = "LCAO_xc_PBE/groundstate"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = groundstate_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'dzp',
                xc='PBE'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
    

    if setting_number == 11:
        name = "LCAO_xc_PBE/excited"
        print(f"{time.ctime()} - Relaxation started for {name} ({setting_number})")
        save_path = 'Assignment 2/problem_8/result/' + name + '/'
        if not os.path.exists(save_path): os.makedirs(save_path)
        run_relaxation(
            xyz_path = excited_path,
            save_path = save_path,
            calc = GPAW(
                nbands = 10,
                h = 0.25,
                txt = save_path + 'gpaw-out.txt',
                occupations = FermiDirac(0.05),
                setups = {'Na': '1'},
                mode = 'lcao',
                basis = 'dzp',
                xc='PBE'
            )
        )
        print(f"{time.ctime()} - Relaxation done for {name} ({setting_number})")
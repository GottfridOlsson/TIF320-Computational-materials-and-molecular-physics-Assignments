# do a new GPAW calculation with the most stable clusters (6, 7, and 8 atoms)
# save the wavefunctions in cube files
# visualized the wavefunctions (plot only the occupied bands. (Which ones are occupied? Count the electrons or look in the log file from gpaw!))
# Please include this script in your report.

### FLOW-CHART ###

# 1. get xyz for stable clusters
# 2. put them into a database
# 3. load it
# 4. set parameters for GPAW
# 5. calculate the electron density
# 6. calculate wavefunction from the electron density
# 7. save the wavefunction as cube file
# 8. try open the cube file in PyMOL to see if it works


# test: https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/wavefunctions/plotting/plot_wave_functions.html#creating-wave-function-cube-files

from ase.io import write, read
from ase.units import Bohr
from gpaw import restart
from ase import Atoms
from gpaw import GPAW

# Declare paths
problem10_path = './Assignment 2/problem_10'
problem5_path = './Assignment 2/problem_5'
xyz_path_start = f'{problem5_path}/structures/natoms_'
xyz_path_end = '_groundstate_gpaw.xyz'

# atoms = read(xyz_path)


for natoms in [6]:

    # Read stable configuration files
    xyz_path = xyz_path_start + str(natoms) + xyz_path_end
    atoms = read(xyz_path)

    # Set GPAW and calculate
    calc = GPAW(nbands=natoms, h=0.2, txt=None, setups={'Na': '1'})
    atoms.calc = calc
    energy = atoms.get_potential_energy()

    # Save wave functions
    #saved_wavefunction_path = problem10_path + '/wavefunctions/' + str(natoms) + '_.gpw'
    #calc.write(saved_wavefunction_path, mode='all')
    #print(f"wrote wavefunctions for natoms={natoms} to path: {saved_wavefunction_path}")

    # Load binary file and get calculator (this might be unnecessary, but kept while testing stuff)
    #atoms, calc = restart(saved_wavefunction_path)

    # Write wavefunctions to cube files
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        wf = calc.get_pseudo_wave_function(band=band)
        fname = f'Na{natoms}_band_{band}.cube'
        print('writing wf', band, 'to file: ' + problem10_path, fname)
        write(problem10_path + '/cubefiles/'+ fname, atoms, data=wf * Bohr**1.5)
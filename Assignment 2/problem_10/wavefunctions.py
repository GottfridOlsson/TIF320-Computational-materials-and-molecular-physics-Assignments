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

problem10_path = './Assignment 2/problem_10/'

from ase.io import write
from ase.units import Bohr
from gpaw import restart

from ase import Atoms
from gpaw import GPAW

if True:
    d = 1.1   # bondlength of hydrogen molecule
    a = 5.0   # sidelength of unit cell
    c = a / 2
    atoms = Atoms('CO',
                positions=[(c - d / 2, c, c),
                            (c + d / 2, c, c)],
                cell=(a, a, a))

    calc = GPAW(nbands=5, h=0.2, txt=None) #h = 0.2
    atoms.calc = calc

    # Start a calculation:
    energy = atoms.get_potential_energy()

    # Save wave functions:
    calc.write(problem10_path + 'CO.gpw', mode='all')



basename = 'CO'

# load binary file and get calculator
atoms, calc = restart(problem10_path + 'CO.gpw')

# loop over all wfs and write their cube files
nbands = calc.get_number_of_bands()
for band in range(nbands):
    wf = calc.get_pseudo_wave_function(band=band)
    fname = f'{basename}_{band}.cube'
    print('writing wf', band, 'to file: ' + problem10_path, fname)
    write(problem10_path + fname, atoms, data=wf * Bohr**1.5)
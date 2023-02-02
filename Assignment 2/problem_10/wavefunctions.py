from ase.io import write, read
from ase.units import Bohr
from gpaw import GPAW

problem10_path = './Assignment 2/problem_10'
problem5_path = './Assignment 2/problem_5'
xyz_path_start = f'{problem5_path}/structures/natoms_'
xyz_path_end = '_groundstate_gpaw.xyz'

for natoms in [6, 7, 8]:
    print(f"Calculating for natoms={natoms}")

    # Read xyz files for the most stable configurations
    xyz_path = xyz_path_start + str(natoms) + xyz_path_end
    atoms = read(xyz_path)

    # Set calculator as GPAW
    calc = GPAW(nbands=5, h=0.2, setups={'Na': '1'}, mode='lcao', basis='dzp', txt=None) # tried nbands=0, but did not work for all natoms ("nbands=0 will give zero empty bands.") #https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html
    atoms.calc = calc
    energy = atoms.get_potential_energy() #it complains if this is not here (do not know why)

    # Get and write wavefunctions to cube files
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        wavefunction = calc.get_pseudo_wave_function(band=band)
        cube_file_path = problem10_path + f'/cubefiles/Na{natoms}_band_number_{band}.cube'
        print(f'Writing wavefunction for Na{natoms} for band number ' + str(band) + f' to file: {cube_file_path}')
        write(cube_file_path, atoms, data=wavefunction * Bohr**1.5)
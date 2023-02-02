from ase.io import write, read
from ase.units import Bohr
from gpaw import GPAW


for natoms in [6, 7, 8]:
    print(f"Calculating wavefunction for Na{natoms}")

    # Read xyz files for the most stable configurations
    xyz_path = f'./Assignment 2/problem_5/structures/natoms_{natoms}_groundstate_gpaw.xyz'
    atoms = read(xyz_path)

    # Set calculator as GPAW
    calc = GPAW(nbands=5, h=0.2, setups={'Na': '1'}, mode='lcao', basis='dzp', txt=None)
    atoms.calc = calc
    energy = atoms.get_potential_energy()

    # Get and write wavefunctions to cube files
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        wavefunction = calc.get_pseudo_wave_function(band=band)
        cube_file_path = f'./Assignment 2/problem_10/cubefiles/Na{natoms}_band_number_{band}.cube'
        print(f'Writing wavefunction for Na{natoms} for band number {band} to file: {cube_file_path}')
        write(cube_file_path, atoms, data=wavefunction * Bohr**1.5)
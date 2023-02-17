# Adapted from:
# https://wiki.fysik.dtu.dk/ase/ase/vibrations/modes.html
# https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.CrystalThermo.get_entropy

from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.build import molecule

molecule_names = ['CO', 'O2']
symmetry_numbers = [1, 2]
spins = [0, 1] # singlet and triplet for CO and O2 respectively

# Temperature and pressure
temperature = 300 # Kelvin
pressure = 1*1e5 # bar*1e5 = Pa

for i in len(molecule_names):
    
    molecule_name = molecule_names[i]
    symmetry_number = symmetry_numbers[i]
    spin = spins[i]

    # Create molecule
    atoms = molecule(molecule_name)
    atoms.calc = EMT()
    dyn = QuasiNewton(atoms)
    dyn.run(fmax=0.01)
    potentialenergy = atoms.get_potential_energy()

    # Set up and run vibrational analysis
    vib = Vibrations(atoms)
    vib.run()


    # Calculate quantities
    potential_energy = atoms.get_potential_energy() # eV
    vibrational_energies = vib.get_energies() # eV
   
    ideal_gas = IdealGasThermo(vib_energies=vibrational_energies, 
                                geometry='linear',
                                potentialenergy=potential_energy,
                                atoms=atoms,
                                symmetrynumber=symmetry_number,
                                spin=spin)
    entropy = ideal_gas.get_entropy(temperature=temperature,pressure=pressure) # eV/K
# Adapted from:
# https://wiki.fysik.dtu.dk/ase/ase/vibrations/modes.html
# https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.CrystalThermo.get_entropy

import numpy as np
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from ase.io import read
from gpaw import GPAW, PW

def print_arrays_to_CSV(path_to_CSV_file, *args, print_message=False):
    """Prints array(s) with corresponding header(s) to a file with comma separated values (CSV)

        Input:
            path_to_csv: the path to where the CSV file should be printed

            *args: array(s) and corresponding header(s) in this format:
                    header_1, array_1, header_2, array_2, ..., header_n, array_n

            print_message: displays a message "Sucessfully printed CSV file (...)" (default False)

        Output:
            A CSV file with utf-8 formatting at path_to_csv, with the array(s) as column(s) and corresponding header(s)

            Lines larger down than the length of array(s) are printed as ',' without an empty space [Plot-Data once complained when there was an empty space]
        
        Warnings:
            ValueError: if the length of args is not even

            ValueError: if the number of array(s) is not equal to the number of header(s)
    """

    if len(args) % 2 != 0:
        raise ValueError("WARNING: the number of arrays + headers is not even. This may cause errors in printing!")
    
    arrays, headers, lines_per_array = [], [], []

    for index, arg in enumerate(args):

        if index % 2 == 0:
            headers.append(arg)
        else:
            arrays.append(arg)
            lines_per_array.append(len(arg))
      

    if len(arrays) != len(headers):
        raise ValueError("WARNING: the number of arrays does not equal the number of headers!")


    with open(path_to_CSV_file, 'w', encoding="utf-8") as CSV_file:
        
        # Print header line
        for header_number, header in enumerate(headers):
            CSV_file.write(str(header))

            # print comma separator between values (CSV) or newline at end of line
            if header_number != len(headers) - 1:
                CSV_file.write(",")
            else:
                CSV_file.write("\n")

        # Print CSV data
        for line in range(max(lines_per_array)):

            for array_number, array in enumerate(arrays):
                # print value of array at line, or if outside length of aray print nothing
                try:
                    CSV_file.write(str(array[line]))
                except:
                    CSV_file.write(str("")) 
                
                # print comma separator between values (CSV) or newline at end of line
                if array_number != len(arrays) - 1:
                    CSV_file.write(",")
                else:
                    CSV_file.write("\n")
    

    if print_message:
        print(f"Successfully printed {len(arrays)} arrays to CSV file at path: '{path_to_CSV_file}'")


molecule_names = ['CO', 'O2']
symmetry_numbers = [1, 2]
spins = [0, 1] # singlet and triplet for CO and O2 respectively

# Temperature and pressure
temperatures = np.linspace(100, 2000, 1901) # Kelvin
pressure_T5 = 1e5 # 1 bar = 10e5 Pa
pressure_T7 = 101325 #Pa, 1 atmosphere of pressure = 101325 Pa


entropies_experimental_units= []


for i in range(len(molecule_names)):
    entropies = []
    molecule_name = molecule_names[i]
    symmetry_number = symmetry_numbers[i]
    spin = spins[i]

    # Read relaxed structure from T4 and do vibrational analysis
    atoms = read(f"output_T4/{molecule_name}_relaxed_molecule_structure.xyz")
    calc = GPAW(xc='PBE',
            mode=PW(450),
            kpts={'gamma': True},
            txt=f"output_T5/{molecule_name}_calc.txt")
    atoms.set_calculator(calc)
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

    for temperature in temperatures:
        entropy = ideal_gas.get_entropy(temperature=temperature,pressure=pressure_T7) # eV/K
        entropies.append(entropy)

        # Convert to J mol^-1 K^-1
        #eV_2_J = 1.602176565e-19
        #per_mole = 6.0221415e23 # Avogadros constant
        #entropies_experimental_units.append(entropy*eV_2_J*per_mole) # J mol^-1 K^-1

    print_arrays_to_CSV(f"output_T5/TIF320_A4_T5_entropy_vs_temperature_at_P={pressure_T7}Pa_{molecule_name}.csv", 
                        "Temperature (K)", temperatures, 
                        "Entropy in ideal gas approximation (eV/K)", entropies)
        #print_arrays_to_CSV(f"output_T5/TIF320_A4_T5_entropy_T={temperature}K_P={pressure:.2e}Pa_CO_O2.csv", 
        #                    "Molecule", molecule_names, 
        #                    "Entropy in ideal gas approximation (eV/K)", entropies, 
        #                    "Entropy in ideal gas approximation (J mol^-1 K^-1)", entropies_experimental_units,
        #                    print_message=True)

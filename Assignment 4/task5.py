# Adapted from:
# https://wiki.fysik.dtu.dk/ase/ase/vibrations/modes.html
# https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#ase.thermochemistry.CrystalThermo.get_entropy

from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from ase.build import molecule
from ase.io import read
from gpaw import GPAW, PW
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton

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
temperature = 300 # Kelvin
pressure = 1e5 # 1 bar = 10e5 Pa

entropies = []
for i in range(len(molecule_names)):
    
    molecule_name = molecule_names[i]
    symmetry_number = symmetry_numbers[i]
    spin = spins[i]

    # Read relaxed molecule from Task 4 and run vibrational analysis
    #molecule = read(f"Assignment 4/output_T4/{molecule_name}_relaxed_molecule_structure.xyz")
    atoms = molecule(molecule_name)
    atoms.calc = EMT()
    dyn = QuasiNewton(atoms)
    dyn.run(fmax=0.01)
    #potentialenergy = atoms.get_potential_energy()

    vib = Vibrations(atoms)
    vib.run()
    #vib_energies = vib.get_energies()

    #calc = GPAW(xc='PBE',
    #            mode=PW(150),
    #            #kpts={'gamma': True},
    #            txt=f"Assignment 4/output_T5/{molecule_name}_GPAW.txt")
    #molecule.set_calculator(calc)
    #vib = Vibrations(molecule)
    #vib.run()

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
    entropies.append(entropy)

    print_arrays_to_CSV("Assignment 4/output_T5/TIF320_A4_T5_entropy_T={temperature}K_P={pressure}Pa_CO_O2.csv", 
                        "Molecule", molecule_names, 
                        "Entropy in ideal gas approximation (eV/K)", entropies, 
                        print_message=True)

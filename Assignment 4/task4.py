from ase.build import molecule
from ase.io import write
from gpaw import GPAW, PW
from ase.optimize import GPMin

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
cell_length = 12 #angstrom
cutoff_energy = 450 #electron volt
output_path_start = "Assignment 4/output_T4/"
potential_energy = []

for molecule_name in molecule_names:

    # Create molecule as atoms object, center, and set calculator
    atoms = molecule(molecule_name, 
                    cell=(cell_length, cell_length, cell_length))
    atoms.center()
    write(f"{output_path_start}{molecule_name}_initial_molecule_structure.png", atoms, rotation='10z,-80x')
    
    calc = GPAW(xc='PBE',
                mode=PW(cutoff_energy),
                kpts={'gamma': True},
                txt=f"{output_path_start}{molecule_name}_calc.txt")

    atoms.set_calculator(calc)
    dyn = GPMin(atoms, 
                trajectory=f"{output_path_start}GPMin_{molecule_name}.traj", 
                logfile=f"{output_path_start}GPMin_{molecule_name}.log")
    dyn.run(fmax=0.01, steps=50) #0.001, 100

    potential_energy.append(atoms.get_potential_energy())
    write(f"{output_path_start}{molecule_name}_relaxed_molecule_structure.xyz", atoms)

    print_arrays_to_CSV(f"{output_path_start}TIF320_A4_T4_potential-energy_of_CO_O2_gas_relaxed.csv", 
                        "Molecule symbol", molecule_names, 
                        "Potential energy (eV)", potential_energy, 
                        print_message=True)
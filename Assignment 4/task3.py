# Construct surface
# Example from https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/structureoptimization/surface/surface.html
# also: https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html

from ase.build import fcc111
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



elements = ['Au', 'Pt', 'Rh'] 
a_from_T1 = [4.177, 3.970, 3.840] # angstrom
output_path_start = 'Assignment 4/output_T3/'
potential_energy_of_surface = []

for i, element in enumerate(elements):
    print(f"i={i}, element={element}")

    # Construct (111)-surface for Au, Pt, Rh
    # 3-layered 3X3 surface cell, 6 angstrom of vacuum in +-z-direction
    surface = fcc111(element, (3, 3, 3), a=a_from_T1[i], vacuum=6.0)
    write(f"{output_path_start}{element}_111-surface_initialized_angled-view.png", surface, rotation='10z,-80x, 5y')
    continue

    # Calculator
    k = (4, 4, 1) #k-sampling according to problem description
    cutoff_energy = 450 # eV
    calc = GPAW(mode=PW(cutoff_energy),
                kpts=k,
                txt=f"{output_path_start}{element}_111-surface_calc.txt")

    # Relax surface (code and settings from Assignment 2, 'relaxation.py')
    surface.set_calculator(calc)
    dyn = GPMin(surface, 
                trajectory=f"{output_path_start}GPMin_{element}_111-surface.traj", 
                logfile=f"{output_path_start}GPMin_{element}_111-surface.log")
    dyn.run(fmax=0.001, steps=100) #0.001, 100
    write(f"{output_path_start}{element}_111-surface_relaxed.png", surface)
    calc.write(f"{output_path_start}{element}_111-surface_relaxed.gpaw")
    
    # Get surface energy (=potential energy?)
    potential_energy_of_surface.append(surface.get_potential_energy())


    print_arrays_to_CSV("Assignment 4/output_T3/TIF320_A4_T3_surface_energy_Pt-Rh.csv", 
                        "Element symbol", elements, 
                        "Potential energy of surface (eV)", potential_energy_of_surface, 
                        print_message=True)
    


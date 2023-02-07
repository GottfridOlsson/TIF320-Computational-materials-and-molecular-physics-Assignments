
def print_arrays_to_CSV(path_to_CSV_file, *args, print_message=False):
    """Prints array(s) with corresponding header(s) to a file with comma separated values (CSV)

        Input:
            path_to_csv: the path to where the CSV file should be printed

            *args: array(s) and corresponding header(s) in this format:
                    header_1, array_1, header_2, array_2, ..., header_n, array_n

            print_message: displays a message "Sucessfully printed CSV file (...)" (default False)

        Output:
            A CSV file with utf-8 formatting at path_to_csv, with the array(s) as column(s) and corresponding header(s)

            Lines larger down than the length of array(s) are printed with an empty space
        
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
                CSV_file.write(", ")
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
                    CSV_file.write(", ")
                else:
                    CSV_file.write("\n")
    

    if print_message:
        print(f"Successfully printed {len(arrays)} arrays to CSV file at path: '{path_to_CSV_file}'")
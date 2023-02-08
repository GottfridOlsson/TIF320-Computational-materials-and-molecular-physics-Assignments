import numpy as np


# Declare paths
our_path   = 'Assignment 3/TIF320_A3_RDF_histogram_our_simulation.csv'
their_path = 'Assignment 3/TIF320_A3_RDF_histogram_their_simulation.csv'
paths = [our_path, their_path]

integral_values = []
for path in paths:

    # Read RDF histogram from file
    RDF_data = np.loadtxt(path, delimiter=',', skiprows=1)
    radial_coordinate = RDF_data[:,0]
    RDF = RDF_data[:,2]

    # For the first coordination number, scale RDF by 4*pi*r^2 and average density
    RDF = RDF*4*np.pi*radial_coordinate*radial_coordinate * (24/8.956**3) # average density: 24 oxygen atoms per unit cell

    # Declare radial value for minimum of RDF (taken from plot/CSV-file of RDF) and get its index
    r_at_first_minimum = 3.2
    r_index_at_first_minimum = 0

    for r in radial_coordinate:
        if r < r_at_first_minimum:
            r_index_at_first_minimum += 1
        else:
            break

    # Numerical integration of RDF up to r_minimum (first solvation shell)
    integral_values.append(np.trapz(RDF[0:r_index_at_first_minimum], radial_coordinate[0:r_index_at_first_minimum]))

# Print result
print(f"Integral value of our simulation:   {integral_values[0]:.3f}")
print(f"Integral value of their simulation: {integral_values[1]:.3f}\n")


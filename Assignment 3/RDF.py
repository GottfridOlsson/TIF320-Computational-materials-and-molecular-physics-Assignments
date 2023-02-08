import numpy as np
from ase.io.trajectory import Trajectory
import util
import matplotlib.pyplot as plt #


# Read trajectories for cluster with Na
trajectories = Trajectory('Assignment 3/logs/nose_hoover_trajectory.traj')

# Create bins for histogram to save RDF
cell_length = trajectories[0].cell[0][0]
r_min, r_max, number_of_bins = 1.5, cell_length/2, 57
histogram = np.array([0]*number_of_bins)
binwidth = (r_max-r_min) / number_of_bins
radial_coordinates = np.array([r_min+binwidth*(i+0.5) for i in range(number_of_bins)])

trajectory_number = 0
for atoms in trajectories:
    trajectory_number += 1
    if trajectory_number % 500 == 0:
        print(f"Calculating histogram for trajectory number {trajectory_number}")

    # Select index for Na and O atoms
    chemical_symbols = np.array(atoms.get_chemical_symbols())
    Na_atom_index    = np.where(chemical_symbols=='Na')[0]
    O_atoms_indexes  = np.where(chemical_symbols=='O')[0]
    
    # Get distance (angstrom) between atoms using the minimum image convention (mic)
    distances = atoms.get_distances(Na_atom_index, O_atoms_indexes, mic=True) 
    
    # Add distances to histogram to produce radial distribution function (RDF)   
    for distance in distances:
        if distance < r_min or distance > r_max:
            continue

        bin_index = int(np.floor( (distance-r_min)/binwidth ))
        histogram[bin_index] += 1


# Scale histogram by average density and by volume of spherical shells and by number of time steps
histogram_scaled = histogram / ( 24 / cell_length**3) # average density is 24 particles (oxygen) in unit cell
histogram_scaled = histogram_scaled / (4*np.pi*radial_coordinates**2*binwidth)
histogram_scaled = histogram_scaled / trajectory_number

# Save histogram to CSV
util.print_arrays_to_CSV("Assignment 3/TIF320_A3_RDF_histogram.csv", 
                            "Radial coordinate (angstrom)", radial_coordinates, 
                            "Histogram of radial distances (counts)", histogram,
                            "Scaled histogram (counts)", histogram_scaled)

plt.plot(radial_coordinates, histogram_scaled, 'o') #
plt.savefig("RDF_histogram.pdf") #
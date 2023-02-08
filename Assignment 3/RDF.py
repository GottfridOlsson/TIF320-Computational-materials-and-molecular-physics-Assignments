import numpy as np
from ase.io.trajectory import Trajectory
import util


# Read trajectories for cluster with Na
trajectories = Trajectory('Assignment 3/logs/nose_hoover_trajectory.traj')

# Create bins for histogram to save RDF
r_min, r_max, binwidth = 0, 10, 0.05
number_of_bins = int((r_max-r_min)/binwidth)
histogram = [0]*number_of_bins
radial_coordinates = [r_min+binwidth*i for i in range(number_of_bins)]

trajectory_number = 0
for atoms in trajectories:
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
        if distance > r_max:
            continue
        bin_index = int(np.floor(distance/binwidth))
        histogram[bin_index] += 1
    
    trajectory_number += 1


# Save histogram to CSV
util.print_arrays_to_CSV("Assignment 3/RDF_histogram.csv", 
                            "Radial coordinate (angstrom)", radial_coordinates, 
                            "Histogram of distances (count)", histogram,)

  

# TODO:
# From these, calculate (or add) the to the RDF
# Loop this forall time steps
# Print RDF to file (radial distance, RDF value)

# skala med 4*pi*r**2
# div med antal punkter i histogrammet

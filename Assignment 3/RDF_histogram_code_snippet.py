# Create histogram to save the RDF to
cell_length = trajectories[0].cell[0][0]
r_min, r_max, number_of_bins = 1.5, cell_length/2, 57
histogram = np.array([0]*number_of_bins)
binwidth = (r_max-r_min)/number_of_bins
radial_coordinates = np.array([r_min+binwidth*(i+0.5) for i in range(number_of_bins)])

# Set counter for number of trajectories (number of time steps) used
trajectory_number = 0

for atoms in trajectories[first_trajectory:last_trajectory]:
    trajectory_number += 1

    # Select index for Na and O atoms
    chemical_symbols = np.array(atoms.get_chemical_symbols())
    Na_atom_index    = np.where(chemical_symbols=='Na')[0]
    O_atoms_indexes  = np.where(chemical_symbols=='O')[0]
        
    # Get distances between atoms using the minimum image convention (mic)
    distances = atoms.get_distances(Na_atom_index, O_atoms_indexes, mic=True) 
    
    # Add distances to histogram   
    for distance in distances:
        if distance < r_min or distance > r_max:
            continue

        bin_index = int(np.floor((distance-r_min)/binwidth))
        histogram[bin_index] += 1

# Scale histogram by volume of spherical shells, by average density, and by number of time steps used
histogram_scaled = histogram/(4*np.pi*radial_coordinates**2*binwidth)
histogram_scaled = histogram_scaled/(24/cell_length**3) # average density is 24 particles (oxygen) per unit cell
histogram_scaled = histogram_scaled/trajectory_number
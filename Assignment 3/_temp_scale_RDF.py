import numpy as np
import util

RDF_data = np.loadtxt("Assignment 3/TIF320_A3_RDF_histogram_their_simulation_without_Na.csv", delimiter=',', skiprows=1)

radial_coordinates = RDF_data[:,0]
histogram = RDF_data[:,1]
RDF = RDF_data[:,2]

histogram_scaled = RDF/24 #divide by # oxygen atoms

util.print_arrays_to_CSV("Assignment 3/TIF320_A3_RDF_histogram_their_simulation_without_Na.csv", 
                            "Radial coordinate (angstrom)", radial_coordinates, 
                            "Histogram of radial distances (counts)", histogram,
                            "Scaled histogram (counts)", histogram_scaled)
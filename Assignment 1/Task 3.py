import numpy as np
import matplotlib.pyplot as plt
from util import *
import hydrogen as hydrogen

def total_probability_of_radial_wavefunction(wavefunction, radial_coordinates, h):
    wavefunction = np.array(wavefunction)
    wavefunction_squared = wavefunction*wavefunction
    return np.trapz(4*np.pi*radial_coordinates*radial_coordinates*wavefunction_squared, radial_coordinates, h)
    # from radial to spherical wavefunction psi(r): --> 4 * pi * r^2 * psi(r)  (integral over space)

def normalize_radial_wavefunction(wavefunction, radial_1D_grid, h):
    integral = total_probability_of_radial_wavefunction(wavefunction, radial_1D_grid, h)
    return wavefunction / np.sqrt(integral)

def create_diagonal_matrix_from_array(array):
    return np.diag(array)

def reciprocal_of_array_handle_division_by_zero(array):
    c = []
    for i in range(len(array)):
        if array[i] == 0:
            c.append(3e99) #'infinity'
        else:
            c.append(1/array[i])
    return c


def divide_arrays_by_each_other(a, b):
    c = []
    for i in range(len(a)):
        if b[i] == 0:
            c.append(1e99)
        else:
            c.append(a[i]/b[i])
    return c


if __name__ == "__main__":

    # define linspace
    N = 1000
    linspace_start, linspace_end = 0, 10
    h = (linspace_end - linspace_start) / (N-1)

    # create linspaces (small letters) and matrices (big letters) 
    r         = create_discretized_1D_space(linspace_start, linspace_end, N, h)
    D2        = create_matrix_D2_finite_difference(N, h)
    R         = create_diagonal_matrix_from_array(r)
    r_inverse = reciprocal_of_array_handle_division_by_zero(r)
    R_inverse = create_diagonal_matrix_from_array(r_inverse)
    r_not_singular = r[1:-1]
    #V_H  = create_diagonal_matrix_from_array()
    #V_x  = create_diagonal_matrix_from_array()
    #V_c  = create_diagonal_matrix_from_array()


    # Solve for hydrogen as a test #
    psi_hydrogen_theoretical = hydrogen.ground_state_wavefunction(r)
    hydrogen_matrix_H = (-0.5*D2 - R_inverse) #Schrödinger equation: H*psi = E*psi, so eigen 

    # solve for eigenvalues and eigenvectors
    hydrogen_eigenvalues, hydrogen_eigenvectors = np.linalg.eigh(hydrogen_matrix_H)
    E_hydrogen = hydrogen_eigenvalues[1] 
    u_hydrogen = -hydrogen_eigenvectors[:,1] #eigh orders eigenvalues, so eigenvector corresponding to smallest eigenvalue is the ground state, also  eigenvectors are in columns

    # wavefunction of hydrogen
    psi_hydrogen = divide_arrays_by_each_other(u_hydrogen, np.sqrt(4*np.pi)*r) #definition of u (=u_hydrogen), equation (34)
    psi_hydrogen_not_singular = psi_hydrogen[1:-1] #get rid of bad first point
    psi_hydrogen_not_singular = normalize_radial_wavefunction(psi_hydrogen_not_singular, r_not_singular, h)
    psi_hydrogen = normalize_radial_wavefunction(psi_hydrogen, r, h)

    # print to verify reasonability
    print("\nHydrogen")
    print(f"Ground state energy: {E_hydrogen:.6f} (theoretically: {hydrogen.ground_state_energy():.3f} (a.u.) )")
    print(f"Wavefunction's total probability (theoretical): {total_probability_of_radial_wavefunction(psi_hydrogen_theoretical, r, h):.6f} (normalized if 1)")
    print(f"Wavefunction's total probability (with singularity): {total_probability_of_radial_wavefunction(psi_hydrogen, r, h):.6f} (normalized if 1)")
    print(f"Wavefunction's total probability (not singularity): {total_probability_of_radial_wavefunction(psi_hydrogen_not_singular, r_not_singular, h):.6f} (normalized if 1)")



    # PLOT #
    plt.plot(r, psi_hydrogen_theoretical,                   color='black', marker='', linestyle='-', label='Theoretical hydrogen wavefunction')
    plt.plot(r_not_singular, psi_hydrogen_not_singular,     color='red', marker='.', linestyle='', label='Calculated hydrogen wavefunction')
    plt.xlabel('Radial distance r (atomic units)')
    plt.ylabel('Wavefunction')
    plt.grid()
    plt.legend()
    plt.show()



    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #



    output_path_wavefunction = f'Assignment 1/output/A1_Task3_hydrogen_wavefunction_N={N}.csv'  
    with open(output_path_wavefunction,'w') as CSV_file:
        CSV_file.write(f"Radial distance r (atomic units), theoretical hydrogen ground state wavefunction (atomic units), calculated hydrogen ground state wavefunction (atomic units) with {N} points\n")
        for line in range(N):
                CSV_file.write(str(r[line]) + ", " + str(psi_hydrogen_theoretical[line]) + ", " + str(psi_hydrogen[line]) + "\n")

    output_path_energy = f'Assignment 1/output/A1_Task3_hydrogen_energy_N={N}.txt'  
    with open(output_path_energy,'w') as file:
        file.write(f"Calculated ground state energy of hydrogen: {E_hydrogen:.8f} (atomic units)\n")
        file.write(f"Theoretical ground state energy of hydrogen: -0.5 (atomic units)\n")
        file.write(f"Number of points in discretized radial coordinate: {N}")
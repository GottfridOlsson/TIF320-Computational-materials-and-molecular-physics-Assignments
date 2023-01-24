import numpy as np
import matplotlib.pyplot as plt
from util import *
import hydrogen as hydrogen

def epsilon_hom_x(n):
    # n is the electron density
    return -(3/4) * (3*n /np.pi)**(1/3)


def r_s_from_n(n):
    # n is the electron density, n = 3 / (4*pi*r_s^3)
    r_s = (3 / (4 * np.pi * n) )**(1/3)
    return r_s


def epsilon_hom_c(r_s):
    # constants from problem description
    A = 0.0311
    B = -0.048
    C = 0.0020
    D = -0.0116
    gamma = -0.1423
    beta_1 = 1.0529
    beta_2 = 0.3334

    if r_s >= 1:
        return gamma / (1 + beta_1*np.sqrt(r_s) + beta_2*r_s)
    if r_s < 1:
        return A*np.log(r_s) + B + C*r_s*np.log(r_s) + D*r_s


def epsilon_hom_xc(n, r_s):
    return epsilon_hom_x(n) + epsilon_hom_c(r_s)

# Potential
def V_xc(n):
    r_s = r_s_from_n(n)
    eps_xc = epsilon_hom_xc(n, r_s)
    return eps_xc + n * np.gradient(eps_xc, n)

# Compute potential from electron density by solving radial poission equation
def V_H(r, n, self_correlation=False):
    
    # Solve poission equation (this is task 2)
    D2 = create_matrix_D2_finite_difference(N, h)
    u_squared_divided_by_r = 4 * np.pi * n * r
    U_0 = np.linalg.solve(D2, -u_squared_divided_by_r) 

    # Compute potential
    U = U_0 + r/np.max(r)
    V_sH = U/r 

    if self_correlation: return 2 * V_sH
    else:                return V_sH

def solve_kohn_sham()

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
    with open(output_path_wavefunction, 'w') as CSV_file:
        CSV_file.write(f"Radial distance r (atomic units), theoretical hydrogen ground state wavefunction (atomic units), calculated hydrogen ground state wavefunction (atomic units) with {N} points\n")
        for line in range(N):
                CSV_file.write(str(r[line]) + ", " + str(psi_hydrogen_theoretical[line]) + ", " + str(psi_hydrogen[line]) + "\n")

    output_path_energy = f'Assignment 1/output/A1_Task3_hydrogen_energy_N={N}.txt'  
    with open(output_path_energy, 'w') as file:
        file.write(f"Calculated ground state energy of hydrogen: {E_hydrogen:.8f} (atomic units)\n")
        file.write(f"Theoretical ground state energy of hydrogen: -0.5 (atomic units)\n")
        file.write(f"Number of points in discretized radial coordinate: {N}")
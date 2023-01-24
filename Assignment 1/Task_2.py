import numpy as np
import matplotlib.pyplot as plt
from util import *
import hydrogen as hydrogen

# Solve Poisson's equation d2U/dr2 = -u2/r,
# using boundary conditions U(0) = 0, U(r_max) = 1
# r must be regular linspace starting at 0
def solve_poisson(r, u):

    # Get linspace properties
    h = r[1]-r[0]
    N = len(r)

    # Solve Poisson's equation with vanishing boundary condition
    D2 = create_matrix_D2_finite_difference(N, h)
    U_0 = np.linalg.solve(D2, -u**2 / r)

    # Fix boundary conditions
    U = U_0 + r/np.max(r)

    return U

if __name__ == "__main__":

    # Define linspace, excluding r=0 to avoid singularity
    N = 1000
    start, end = 0, 10
    r = np.linspace(start, end, N+1)[1:]
    h = r[1] - r[0]

    # Get theoretical ground state electron density for hydrogen
    electron_density = hydrogen.ground_state_electron_density(r)

    # Solve Poisson's equation to obtain the potential
    u = np.sqrt(4 * np.pi * electron_density) * r
    U = solve_poisson(r, u)
    V_sH = U/r 

    # Get theoretical result for comparison
    V_Hartree = hydrogen.hartree_potential(r)

    # PLOT #
    plt.plot(r, V_Hartree,  color='black', marker='', linestyle='-', label='Theoretical Hartree potential')
    plt.plot(r, V_sH,       color='red', marker='.', linestyle='', label='Calculated Hartree potential')

    plt.xlabel('Radial distance r (atomic units)')
    plt.ylabel('Energy of potential (atomic units)')
    plt.grid()
    plt.legend()
    plt.show()


    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
    print_arrays_to_CSV(f'Assignment 1/output/A1_Task2_Hartree_potential_N={N}.csv', 
                        "Radial distance r (atomic units)", r, 
                        "theoretical Hartree potential (atomic units)", V_Hartree, 
                        f"calculated Hartree potential (atomic units) with {N} points", V_sH, 
                        print_message=True)
import numpy as np
import matplotlib.pyplot as plt
from util import *
import hydrogen as hydrogen


if __name__ == "__main__":

    N = 100
    linspace_start, linspace_end = 0, 10
    h = (linspace_end - linspace_start)/(N-1)
    r = create_discretized_1D_space(linspace_start, linspace_end, N, h)

    V_Hartree = hydrogen.hartree_potential(r)
    D2 = create_matrix_D2_finite_difference(N, h)
    electron_density = hydrogen.ground_state_electron_density(r)
    u_squared_divided_by_r = 4 * np.pi * electron_density * r

    #solve matrix equation: A*U_0 = b
    U_0 = np.linalg.solve(D2, -u_squared_divided_by_r) 
    # where A is the matrix D2 defined in code (finite difference of second derivative)
    # U_0 is the potential U_0 with boundary conditions U_0(0) = U_0(r_max) = 0,
    #     and U(r) = U_0(r) + r/r_max
    # b is - u**2/r where h is the step in the discretization of grid and
    # u is np.sqrt(4*np.pi*n(r))*r, and n(r) is the electron density for hydrogen (one electron)

    # Get the sought potential from definition of U_0 and relation between V_sH and U
    U = U_0 + r/np.max(r)
    V_sH = U/r 


    # PLOT #
    plt.plot(r, V_Hartree,  color='black', marker='', linestyle='-', label='Theoretical Hartree potential')
    plt.plot(r, V_sH,       color='red', marker='.', linestyle='', label='Calculated Hartree potential')

    plt.xlabel('Radial distance r (atomic units)')
    plt.ylabel('Energy of potential (atomic units)')
    plt.grid()
    plt.legend()
    plt.show()


    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
    output_path = f'Assignment 1/output/A1_Task2_Hartree_potential_N={N}.csv'
    with open(output_path,'w') as CSV_file:
        CSV_file.write(f"Radial distance r (atomic units), theoretical Hartree potential (atomic units), calculated Hartree potential (atomic units) with {N} points\n")
        for line in range(N):
                CSV_file.write(str(r[line]) + ", " + str(V_Hartree[line]) + ", " + str(V_sH[line]) + "\n")
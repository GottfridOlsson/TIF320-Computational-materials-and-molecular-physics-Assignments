import numpy as np
import matplotlib.pyplot as plt
from util import *
import hydrogen as hydrogen

# Solves radial Kohn-Sham equations for a given potential,
# (- 0.5 * d2/dr2 + potential) u = eps u
# returns eps and u for the lowest energy solution.
# u will be normalized to have unity integral of u^2
def solve_kohn_sham(r, potential):
    
    h = r[1]-r[0]
    N = len(r)
    
    D2 = create_matrix_D2_finite_difference(N, h)
    potential_matrix = np.diag(potential)

    # Solve Kohn-Sham matrix equation
    eps_vec, u_mat = np.linalg.eigh(-0.5 * D2 + potential_matrix)

    # eigh will sort eigenvalues in ascending order
    eps = eps_vec[0] 
    u = u_mat[:,0]

    if u[0] < 0: u = -u
    norm2 = np.trapz(u**2, r)
    u *= 1 / np.sqrt(norm2)

    return eps, u


if __name__ == "__main__":

    # Define linspace, excluding r=0 to avoid singularity
    N = 100
    linspace_start, linspace_end = 0, 10
    r = np.linspace(linspace_start, linspace_end, N+1)[1:]
    h = r[1]-r[0]

    # Solve for hydrogen as a test #
    potential = - 1 / r
    E_hydrogen, u_hydrogen = solve_kohn_sham(r, potential)

    # Compute wavefunction from u
    psi_hydrogen = u_hydrogen / (np.sqrt(4*np.pi) * r)
    psi_hydrogen = normalize_radial_wavefunction(psi_hydrogen, r, h)

    # Compare with theoretical
    psi_hydrogen_theoretical = hydrogen.ground_state_wavefunction(r)

    # print to verify reasonability
    print("\nHydrogen")
    print(f"Ground state energy: {E_hydrogen:.6f} (theoretically: {hydrogen.ground_state_energy():.3f} (a.u.) )")
    print(f"Wavefunction's total probability (theoretical): {total_probability_of_radial_wavefunction(psi_hydrogen_theoretical, r, h):.6f} (normalized if 1)")
    print(f"Wavefunction's total probability (calculated): {total_probability_of_radial_wavefunction(psi_hydrogen, r, h):.6f} (normalized if 1)")


    # PLOT #
    plt.plot(r, psi_hydrogen_theoretical, color='black', marker='', linestyle='-', label='Theoretical hydrogen wavefunction')
    plt.plot(r, psi_hydrogen, color='red', marker='.', linestyle='', label='Calculated hydrogen wavefunction')
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
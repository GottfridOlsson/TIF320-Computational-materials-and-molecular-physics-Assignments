import numpy as np
import matplotlib.pyplot as plt
from util import *
from Task_2 import solve_poisson
from Task_3 import solve_kohn_sham

# Compute ground state energy for helium in local density approximation
def compute_energy(r, u, eps, V_H=None, V_xc=None, eps_xc=None):

    local_energy = np.zeros_like(r)
    if not V_H is None:     local_energy += 0.5 * V_H
    if not V_xc is None:    local_energy += V_xc
    if not eps_xc is None:  local_energy -= eps_xc

    E0 = 2 * eps - 2 * np.trapz(local_energy * u**2, r)
    return E0

if __name__ == "__main__":

    # Define linspace, excluding r=0 to avoid singularity
    N = 1000
    linspace_start, linspace_end = 0, 10
    r = np.linspace(linspace_start, linspace_end, N+1)[1:]
    h = r[1]-r[0]

    # Make initial guess for density
    n_s_initial = np.ones_like(r) / (linspace_end - linspace_start)
    u = np.sqrt(4 * np.pi * n_s_initial) * r

    for i in range(0, 100):

        # Compute Hartee potential for current density
        U = solve_poisson(r, u)
        V_sH = U / r
        V_H = V_sH # No self interaction

        # Solve Kohn sham for this potential to get updated density
        potential = -2.0/r + V_H
        eps, u = solve_kohn_sham(r, potential)

        E0 = compute_energy(r, u, eps, V_H=V_H)
        print(E0)

    # Convert u to wavefunction
    psi = u / (np.sqrt(4 * np.pi) * r)

    # print to verify reasonability
    print("\nHydrogen")
    print(f"Ground state energy: {E0:.6f} (a.u.)")
    print(f"Wavefunction's total probability: {total_probability_of_radial_wavefunction(psi, r, h):.6f} (normalized if 1)")
    
    # PLOT #
    plt.plot(r, psi, color='black', marker='', linestyle='-', label='Computed helium wavefunction')
    plt.xlabel('Radial distance r (atomic units)')
    plt.ylabel('Wavefunction')
    plt.grid()
    plt.legend()
    plt.show()


    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
    print_arrays_to_CSV(f'Assignment 1/output/A1_Task4_helium_wavefunction_N={N}.csv', 
                        "Radial distance r (a.u.)", r, 
                        "Calculated helium wavefunction", psi,  
                        print_message=True)
    output_path_wavefunction = f'Assignment 1/output/A1_Task4_helium_wavefunction_N={N}.csv'  

    output_path_energy = f'Assignment 1/output/A1_Task4_helium_energy_N={N}.txt'  
    with open(output_path_energy, 'w') as file:
        file.write(f"Calculated ground state energy of helium: {E0:.8f} (a.u.)\n")
        file.write(f"Number of points in discretized radial coordinate: {N}")
import numpy as np
import matplotlib.pyplot as plt
from util import *
from Task_2 import solve_poisson
from Task_3 import solve_kohn_sham
from Task_4 import compute_energy


def n_s_from_u(r, u):
    return u**2 / (4 * np.pi * r**2)

def epsilon_hom_x(n):
    # n is the electron density
    return -(3/4) * (3 * n / np.pi)**(1/3)

def epsilon_hom_c(n):
    # constants from problem description
    A = 0.0311
    B = -0.048
    C = 0.0020
    D = -0.0116
    gamma = -0.1423
    beta_1 = 1.0529
    beta_2 = 0.3334

    r_s = (3 / (4 * np.pi * n) )**(1/3)
    return np.where(r_s >=1, 
        gamma / (1 + beta_1*np.sqrt(r_s) + beta_2*r_s),
        A*np.log(r_s) + B + C*r_s*np.log(r_s) + D*r_s
    )

def epsilon_hom_xc(n):
    return epsilon_hom_x(n) + epsilon_hom_c(n)

if __name__ == "__main__":

    # Define linspace, excluding r=0 to avoid singularity
    N = 1000
    linspace_start, linspace_end = 0, 10
    r = np.linspace(linspace_start, linspace_end, N+1)[1:]
    h = r[1]-r[0]

    # Make initial guess for density
    n_s_initial = np.exp(-r)
    u = np.sqrt(4 * np.pi * n_s_initial) * r

    for i in range(0, 100):

        # Compute Hartee potential for current density
        U = solve_poisson(r, u)
        V_sH = U / r
        V_H = 2 * V_sH # Include self interaction

        # Compute electron density
        n = n_s_from_u(r, u)

        # Compute exchange correlation potential
        eps_xc = 2 * epsilon_hom_xc(n)
        V_xc = eps_xc + n * np.gradient(eps_xc, n)

        # Solve Kohn sham for this potential to get updated density
        potential = -2.0/r + V_H + V_xc
        eps, u = solve_kohn_sham(r, potential)

        E0 = compute_energy(r, u, eps, V_H=V_H, V_xc=V_xc, eps_xc=eps_xc)
        print(f"Iteration {i+1}, E_0 = {E0:.6f} (a.u.)")

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
    print_arrays_to_CSV(f'Assignment 1/output/A1_Task5_helium_wavefunction_N={N}.csv', 
                        "Radial distance r (a.u.)", r, 
                        "Calculated helium wavefunction", psi,  
                        print_message=True)

    with open(f'Assignment 1/output/A1_Task4_helium_energy_N={N}.txt', 'w') as file:
        file.write(f"Calculated ground state energy of helium: {E0:.8f} (a.u.)\n")
        file.write(f"Number of points in discretized radial coordinate: {N}")
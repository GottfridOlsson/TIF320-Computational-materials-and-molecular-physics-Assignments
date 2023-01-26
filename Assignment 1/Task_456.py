import numpy as np
import matplotlib.pyplot as plt
from util import *
from Task_2 import solve_poisson
from Task_3 import solve_kohn_sham

def wavefunction_anzats(r):
    alpha  = np.array([0.297104, 1.236745, 5.749982, 38.216677])
    C      = np.array([-0.14687613, -0.39315207, -0.41119776, -0.26200681])

    psi = np.zeros_like(r)

    for p in range(0, len(C)):
        psi += C[p] * np.exp( - alpha[p] * r**2 )
    
    return psi

def n_s_from_u(r, u):
    return u**2 / (4 * np.pi * r**2)

def eps_exchange(n):
    return - (3/4) * (3 * n / np.pi)**(1/3)

def V_exchange(n):
    return - (3 * n / np.pi)**(1/3)

def eps_correlation(n):
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

def V_correlation(n):
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
        gamma * 
        (3 + 3.5 * beta_1 * np.sqrt(r_s) + 4 * beta_2 * r_s) /
        (3 * (1 + beta_1 * np.sqrt(r_s) + beta_2 * r_s))
        ,
        A * (np.log(r_s) - 1/3) + 
        B + 
        C * (2 * r_s * np.log(r_s) - r_s) / 3 + 
        2 * D * r_s / 3
    )

# Solves for self consistent wavefunction for one specific rmax and grid size
def find_scf_wavefunction(
    rmax, N, 
    include_exchange = False, 
    include_correlation = False
):
    
    # Create linspace
    r = np.linspace(0, rmax, N+1)[1:]

    # Make initial guess for density
    u = np.sqrt(4 * np.pi) * r * wavefunction_anzats(r)

    # Iterate until convergence to find self consistent solution
    E = 0
    E_old = 0
    for i in range(0, 100):
        
        # Compute Hartee potential for current density
        U = solve_poisson(r, u)
        V_sH = U / r
        
        # If exchange interaction is not included, remove self interaction
        if (include_exchange): V_H = 2 * V_sH 
        else:                  V_H = V_sH

        # Compute electron density
        n = 2 * n_s_from_u(r, u)

        # Compute exchange potential
        eps_xc = np.zeros_like(r)
        V_xc = np.zeros_like(r)
        if (include_exchange): 
            eps_xc += eps_exchange(n)
            V_xc += V_exchange(n)
        if (include_correlation): 
            eps_xc += eps_correlation(n)
            V_xc += V_correlation(n)

        # Solve Kohn sham for this potential to get updated density
        potential = - 2.0 / r + V_H + V_xc
        eps, u = solve_kohn_sham(r, potential)

        # Compute ground state energy
        E_old = E
        E = 2 * eps - 2 * np.trapz(u**2 * (
            0.5 * V_H + V_xc - eps_xc
        ), r)

        print(f"Iteration {i+1}, E_0 = {E:.6f} (a.u.)")

        # Break if energy step is less than 10^-5 eV
        if (27.211 * abs(E - E_old) < 1e-5): break

    # Convert u to wavefunction
    psi = u / (np.sqrt(4 * np.pi) * r)

    return E, r, psi


if __name__ == "__main__":

    test_rmax_convergence = False
    test_N_convergence = False

    include_exchange = True
    include_correlation = True

    if (test_rmax_convergence):
        # Increase r_max until convergence
        rmax_range = np.array([5, 10, 15, 20, 25, 30], dtype=np.float)
        E_range = np.zeros_like(rmax_range)
        grid_density = 100 # points per unit r
        for i, rmax in enumerate(rmax_range):

            N = int(rmax * grid_density)
            print(f"\nr_max = {rmax}, N={N}")
            E, r, psi = find_scf_wavefunction(
                rmax, N,
                include_exchange=include_exchange,
                include_correlation=include_correlation
            )
            E_range[i] = E

        #plt.plot(rmax_range, E_range, "x")
        #plt.xlabel("Grid size $r_{max}$")
        #plt.ylabel("Ground state energy $E0$")
        #plt.grid()
        #plt.show()

        print_arrays_to_CSV(f'Assignment 1/output/A1_Task4_rmax_convergence_density={grid_density}.csv', 
                            "Maximum radial distance r_max (a.u.)", rmax_range, 
                            "Calculated ground state energy", E_range,  
                            print_message=True)

    # Increase N until convergence
    if (test_N_convergence):
        N_range = np.array([750, 1500, 3000, 4500, 6000])
        E_range = np.zeros_like(N_range, dtype=np.float)
        rmax = 30
        for i, N in enumerate(N_range):
            print(f"\nN = {N}")
            E, r, psi = find_scf_wavefunction(
                rmax, N, 
                include_exchange=include_exchange,
                include_correlation=include_correlation
            )
            E_range[i] = E

        #plt.plot(N_range, E_range, "x")
        #plt.xlabel("Number of grid points $N$")
        #plt.ylabel("Ground state energy $E0$")
        #plt.ion()
        #plt.show(block=False)

        print_arrays_to_CSV(f'Assignment 1/output/A1_Task4_N_convergence_rmax={rmax}.csv', 
                            "Number of grid points N", N_range, 
                            "Calculated ground state energy", E_range,  
                            print_message=True)


    # Final calculation with converged values
    rmax = 30
    N = 6000
    print("Computing ground state")
    E, r, psi = find_scf_wavefunction(
        rmax, N, 
        include_exchange=include_exchange,
        include_correlation=include_correlation
    )


    # print to verify reasonability
    print("\nHelium")
    print(f"Ground state energy: {E:.6f} (a.u.)")
    print(f"Wavefunction's total probability: {total_probability_of_radial_wavefunction(psi, r):.6f} (normalized if 1)")
    
    # PLOT #
    plt.plot(r, psi, color='black', marker='', linestyle='-', label='Computed helium wavefunction')
    plt.xlabel('Radial distance r (atomic units)')
    plt.ylabel('Wavefunction')
    plt.grid()
    plt.legend()
    plt.show()


    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
    task = "X"
    if (not include_correlation and not include_exchange): task = "4"
    elif (not include_correlation and include_exchange): task = "5"
    elif (include_correlation and include_exchange): task = "6"
    print_arrays_to_CSV(f'Assignment 1/output/A1_Task{task}_helium_wavefunction_N={N}_rmax={rmax}_x={include_exchange}_c={include_correlation}.csv', 
                        "Radial distance r (a.u.)", r, 
                        "Calculated helium wavefunction", psi,  
                        print_message=True)

    with open(f'Assignment 1/output/A1_Task{task}_helium_energy_N={N}_rmax={rmax}_x={include_exchange}_c={include_correlation}.txt', 'w') as file:
        file.write(f"Calculated ground state energy of helium: {E:.8f} (a.u.)\n")
        file.write(f"Number of points in discretized radial coordinate: {N}")
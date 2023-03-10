import numpy as np
import matplotlib.pyplot as plt
from util import *

# Wavefunction anzats used in this task: sum of four gaussians
def wavefunction_anzats(r, C, alpha):    
    result = np.zeros_like(r)

    for p in range(0, len(C)):
        result += C[p] * np.exp( - alpha[p] * r**2 )
    
    return result


def normalize_coefficients(C, S):
    n = len(C)
    norm2 = 0
    for p in range(0, n):
        for q in range(0, n):
            norm2 += C[p] * S[p][q] * C[q]
    
    return C / np.sqrt(norm2)


def create_h_matrix(alpha):
    n = len(alpha)
    h = np.zeros((n, n))

    for p in range(0, n):
        for q in range(0, n):

            alpha_sum = alpha[p]+alpha[q]
            h[p][q] = \
                3 * np.pi * alpha[q]    * np.sqrt(np.pi / alpha_sum**3) - \
                3 * np.pi * alpha[q]**2 * np.sqrt(np.pi / alpha_sum**5) - \
                4 * np.pi / alpha_sum
    
    return h
                

def create_S_matrix(alpha):
    n = len(alpha)
    S = np.zeros((n, n))

    for p in range(0, n):
        for q in range(0, n):

            S[p][q] = (np.pi / (alpha[p] + alpha[q]))**1.5

    return S


def create_Q_tensor(alpha):
    n = len(alpha)
    Q = np.zeros((n,n,n,n))

    for p in range(0, n):
        for r in range(0, n):
            for q in range(0, n):
                for s in range(0, n):

                    Q[p][r][q][s] = 2 * np.pi**2.5 / (
                        (alpha[p] + alpha[q]) * (alpha[r] + alpha[s]) *
                        np.sqrt(alpha[p] + alpha[q] + alpha[r] + alpha[s])
                    )

    return Q


def create_F_matrix(C, h, Q):
    n = len(C)
    F = np.zeros((n, n))

    for p in range(0, n):
        for q in range(0, n):

            F[p][q] = h[p][q]

            for r in range(0, n):
                for s in range(0, n):

                    F[p][q] += Q[p][r][q][s] * C[r] * C[s]

    return F


def compute_energy(C, h, Q):
    n = len(C)
    E = 0
    for p in range(0, n):
        for q in range(0, n):

            E += 2 * C[p] * C[q] * h[p][q]

            for r in range(0, n):
                for s in range(0, n):

                    E += Q[p][r][q][s] * C[p] * C[q] * C[r] * C[s]

    return E


if __name__ == "__main__":
    alpha  = np.array([0.297104, 1.236745, 5.749982, 38.216677])

    h = create_h_matrix(alpha)
    S = create_S_matrix(alpha)
    Q = create_Q_tensor(alpha)


    # Pick initial value:
    C = np.array([1, 1, 1, 1])
    C = normalize_coefficients(C, S)
    E = compute_energy(C, h, Q)
    old_E = E

    # Iterate to find self-consistent solution
    for i in range(0, 1000):
        F = create_F_matrix(C, h, Q)
        Eprime, C = solve_generalized_eigenvalue_problem(F, S)
        C = normalize_coefficients(C, S)
        
        old_E = E
        E = compute_energy(C, h, Q)

        # Break if energy step is less than 10^-5 eV
        if (27.211 * abs(E - old_E) < 1e-5): break


    print(f"Ground state energy: {E:.7f} (should be -2.8551716)")
    print(f"C-parameters: {C}")

    N = 1000
    r_lin = np.linspace(0, 5, N)
    phi = np.abs(wavefunction_anzats(r_lin, C, alpha))
   
    # PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
    print_arrays_to_CSV(f'Assignment 1/output/A1_Task1_helium_wavefunction_4gaussians_as_basis_N={N}.csv', 
                        "Radial distance r (atomic units)", r_lin, 
                        f"Calculated helium wavefunction from basis of four gaussians with {N} points (probability amplitude)", phi,
                        print_message=True)

    with open('Assignment 1/output/A1_Task1_helium_wavefunction_alpha_and_C-parameters.txt','w') as file:
        file.write(f"Calculated C-parameters from ground state wavefunction of helium with four gassians as a basis, each on the form: C * exp( - alpha * r^2 )\n")
        file.write(f"alpha[0]: {alpha[0]}\n")
        file.write(f"alpha[1]: {alpha[1]}\n")
        file.write(f"alpha[2]: {alpha[2]}\n")
        file.write(f"alpha[3]: {alpha[3]}\n")
        file.write(f"C[0]: {C[0]}\n")
        file.write(f"C[1]: {C[1]}\n")
        file.write(f"C[2]: {C[2]}\n")
        file.write(f"C[3]: {C[3]}\n")
        file.write(f"Number of points in discretized radial coordinate: {N}\n")
        file.write(f"Ground state energy of helium: {E:.7f} (a.u.)")

    plt.plot(r_lin, phi)
    plt.xlabel("Radial distance from nucleus (a.u.)")
    plt.ylabel("Wavefunction")
    plt.show()
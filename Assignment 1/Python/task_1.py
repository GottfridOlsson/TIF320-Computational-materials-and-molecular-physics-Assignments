import numpy as np
import matplotlib.pyplot as plt

def wavefunction_anzats(r, C, alpha):    
    result = np.zeros_like(r)

    for p in range(0, len(C)):
        result += C[p] * np.exp( - alpha[p] * r**2 )
    
    return result


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


def normalize_coefficients(C, S):
    
    n = len(C)
    norm2 = 0
    for p in range(0, n):
        for q in range(0, n):
            norm2 += C[p] * S[p][q] * C[q]
    
    return C / np.sqrt(norm2)


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


# Solve generalized eigenvalue problem:
#   (F - Eprime S)C = 0   <=>   F C = Eprime S C
# where F and S are matrices, C is an eigenvector and E is an eigenvalue.
# S is assumed to be a real symmetric matrix (diagonalizable)

# Returns smallest eigenvalue and corresponding eigenvector.
def solve_generalized_eigenvalue_problem(F, S):
    
    # Diagonalize S
    d, U = np.linalg.eigh(S)

    # V will have the property V.T @ S @ V = I
    V = U @ np.diag(d ** -0.5)

    # It can be shown that if Cprime is eigenvalue to V.T @ F @ V,
    # then C is V @ Cprime
    Eprime, Cprime = np.linalg.eig(V.T @ (F @ V))
    C = V @ Cprime

    # np.eig does not give ordered eigenvalues! (but np.eigh does)
    i  = np.argmin(Eprime) # index of smallest eigenvalue

    return Eprime[i], C[:, i]

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


r_lin = np.linspace(0, 5, 1000)
phi = np.abs(wavefunction_anzats(r_lin, C, alpha))
plt.plot(r_lin, phi)
plt.xlabel("Radial distance from nucleus (a.u.)")
plt.ylabel("Wavefunction")
plt.show()
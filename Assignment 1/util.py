import numpy as np

# Solve generalized eigenvalue problem:
# (F - Eprime S)C = 0   <=>   F C = Eprime S C
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


# Create linspace for finite difference method
def create_discretized_1D_space(start, end, number_of_points, distance_between_points=None):
    if distance_between_points is not None: 
        h = distance_between_points
    else:
        h = (end-start)/(number_of_points-1)

    linspace = np.array([start + i*h for i in range(number_of_points)])

    return linspace

# Create matrix representation of second derivative operator
def create_matrix_D2_finite_difference(number_of_points_in_discretized_1D_grid, h):
    N = number_of_points_in_discretized_1D_grid

    D2 = np.zeros((N,N))
    i,j = np.indices(D2.shape)

    # Operator matrix for numerical second derivative
    # formula: d2y/dx2 = (y(k-1) - 2y(k) + y(k+1)) / dx2
    D2[i==j]        = -2 / h**2    
    D2[abs(i-j)==1] =  1 / h**2

    return D2

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
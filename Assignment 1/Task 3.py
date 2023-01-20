import numpy as np
import matplotlib.pyplot as plt
import os

def get_current_absolute_path():
    return os.path.abspath(os.path.dirname(__file__))


def create_discretized_1D_space(start, end, number_of_points, distance_between_points=None):
    if distance_between_points is not None: 
        h = distance_between_points
    else:
        h = (end-start)/(number_of_points-1)

    linspace = [start + i*h for i in range(number_of_points)]

    return linspace


def create_matrix_D2_finite_difference(number_of_points_in_discretized_1D_grid, h):
    N = number_of_points_in_discretized_1D_grid

    D2 = np.zeros((N,N))
    i,j = np.indices(D2.shape)

    # Operator matrix for numerical second derivative
    # formula: d2y/dx2 = (y(k-1) - 2y(k) + y(k+1)) / dx2
    D2[i==j]        = -2 / h**2    
    D2[abs(i-j)==1] =  1 / h**2

    return D2

def create_diagonal_matrix_from_array(array):
    #element matrix[i][i] is set to array[i]
    reversed_array = np.flip(array)
    return np.diag(reversed_array)


def normalized_ground_state_wavefunction_hydrogen(radial_distance):
    # in atomic units, Bohr radius a0 = 1
    a0 = 1

    r = np.array(radial_distance)
    f = (1 / np.sqrt(np.pi)) * (1 / a0**(3/2)) * np.exp(-r / a0)
    return f

def ground_state_energy_hydrogen():
    # in atomic units
    return -(1/2)



## MAIN ##

N = 100
linspace_end, linspace_start = 0, 10
h = (linspace_end - linspace_start) / (N-1)

r    = create_discretized_1D_space(linspace_start, linspace_end, N, h)
D2   = create_matrix_D2_finite_difference(N, h)
R    = create_diagonal_matrix_from_array(r)
r_inverse = np.array([1/x for x in r[0:-1]])
r_inverse = np.insert(r_inverse, 0, 10e10, axis=0) #to avoid division by zero
R_inverse = create_diagonal_matrix_from_array(r_inverse)
#V_H  = create_diagonal_matrix_from_array()
#V_x  = create_diagonal_matrix_from_array()
#V_c  = create_diagonal_matrix_from_array()


# Solve for hydrogen as a test #
psi_hydrogen_theoretical = normalized_ground_state_wavefunction_hydrogen(r)
hydrogen_matrix_H = (-0.5*D2 - R_inverse) #Schrödinger equation: H*psi = E*psi, so eigen 


hydrogen_eigenvalues, hydrogen_eigenvectors = np.linalg.eigh(hydrogen_matrix_H)
print(hydrogen_eigenvalues)
E_hydrogen = hydrogen_eigenvalues[1] 
u_hydrogen = hydrogen_eigenvectors[:,1] #eig orderes eigenvalues, so eigenvector corresponding to smallest eigenvalue is the ground state, also  
print(u_hydrogen)
psi_hydrogen = [u / (np.sqrt(4 * np.pi) * r_i) for (u, r_i) in zip(u_hydrogen, r)] #definition of u, equation (34)


print("Hydrogen")
print(f"Ground state energy: {E_hydrogen} (theoretically: {ground_state_energy_hydrogen()} (a.u.) )")




# PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
current_absolute_path = get_current_absolute_path()
output_string = f'\output\A1_Task3_hydrogen_N={N}.csv'
with open(current_absolute_path+output_string,'w') as CSV_file:
    CSV_file.write(f"Radial distance r (atomic units), theoretical hydrogen ground state wavefunction (atomic units), calculated hydrogen ground state wavefunction (atomic units) with {N} points\n")
    #for line in range(N):
    #        CSV_file.write(str(r[line]) + ", " + str(V_Hartree[line]) + ", " + str(V_sH[line]) + "\n")
        



# PLOT #
plt.plot(r, psi_hydrogen_theoretical,  color='black', marker='', linestyle='-', label='Theoretical hydrogen wavefunction')
plt.plot(r, psi_hydrogen,       color='red', marker='.', linestyle='', label='Calculated hydrogen wavefunction')

plt.xlabel('Radial distance r (atomic units)')
plt.ylabel('Wavefunction')
plt.grid()
plt.legend()
plt.show()

# EOF # 
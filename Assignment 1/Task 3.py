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

def normalized_ground_state_wavefunction_hydrogen():
    # in atomic units, Bohr radius a0 = 1
    a0 = 1

    f = 
    return f

## MAIN ##

N = 10
linspace_end, linspace_start = 0, 10
h = (linspace_end - linspace_start) / (N-1)
r = create_discretized_1D_space(linspace_start, linspace_end, N, h)



eigenvalues, eigenvectors = np.linalg.eig()



# PRINT DATA TO CSV FOR PLOT IN PLOT-DATA #
current_absolute_path = get_current_absolute_path()
output_string = f'\output\A1_Task2_N={N}.csv'
with open(current_absolute_path+output_string,'w') as CSV_file:
    CSV_file.write(f"Radial distance r (atomic units), theoretical Hartree potential (atomic units), calculated Hartree potential (atomic units) with {N} points\n")
    for line in range(N):
            CSV_file.write(str(r[line]) + ", " + str(V_Hartree[line]) + ", " + str(V_sH[line]) + "\n")



# PLOT #
plt.plot(r, V_Hartree,  color='black', marker='', linestyle='-', label='Theoretical Hartree potential')
plt.plot(r, V_sH,       color='red', marker='.', linestyle='', label='Calculated Hartree potential')

plt.xlabel('Radial distance r (atomic units)')
plt.ylabel('Energy of potential (atomic units)')
plt.grid()
plt.legend()
plt.show()

# EOF # 
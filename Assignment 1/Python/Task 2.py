import numpy as np
import matplotlib.pyplot as plt



def create_discretized_1D_space(start, end, number_of_points, distance_between_points=None):
    if distance_between_points is not None: 
        h = distance_between_points
    else:
        h = (end-start)/(number_of_points-1)

    linspace = [start + i*h for i in range(number_of_points)]

    return linspace


def Hartree_potential(radial_distance):
    r = np.array(radial_distance)
    return (1.0/r) - (1.0 + 1.0/r)*np.exp(-2.0*r)


def create_matrix_D2_task_2(number_of_points_in_discretized_1D_grid):
    N = number_of_points_in_discretized_1D_grid

    D2 = np.zeros((N,N))
    i,j = np.indices(D2.shape)

    # Operator matrix for numerical second derivative
    # formula: d2y/dx2 = (y(k-1) - 2y(k) + y(k+1)) / dx2
    D2[i==j-1] = 1 / h**2
    D2[i==j]   = -2 / h**2
    D2[i==j+1] = 1 / h**2

    return D2


def ground_state_electron_density_for_hydrogen(radial_distance):
    # atomic units, Bohr radius: a0 = 1
    r  = np.array(radial_distance)
    a0 = 1
    n  = np.exp(-2*r/a0) / (np.pi * a0**3)
    # source: https://en.wikipedia.org/wiki/Hydrogen_atom
    return n





## MAIN ##

linspace_start, linspace_end = 0, 10
N = 1000
h = (linspace_end - linspace_start)/(N-1)
r = create_discretized_1D_space(linspace_start, linspace_end, N, h)

V_Hartree = Hartree_potential(r)
D2 = create_matrix_D2_task_2(N)
electron_density = ground_state_electron_density_for_hydrogen(r)
#u = np.sqrt(4*np.pi*electron_density)*r
#u_squared = 4*np.pi*electron_density*r**2
u_squared_divided_by_r = 4 * np.pi * electron_density * r

U_0 = np.linalg.solve(D2, -u_squared_divided_by_r) #solve matrix equation: A*U_0 = b
# where A is the matrix D2 defined in code (finite difference of second derivative)
# U_0 is the potential U_0 with boundary conditions U_0(0) = U_0(r_max) = 0,
#     and U(r) = U_0(r) + r/r_max is the sought potential (I think?)
# b is - u**2/r where h is the step in the discretization of grid and
# u is np.sqrt(4*np.pi*n(r))*r, and n(r) is the electron density for hydrogen (one electron)

U = U_0 + r/np.max(r) #definition of U_0
# Scale the solution as to satisfy boundary conditions:
#U[0] = 0, U[-1] = 1
#U += 10 #to shift it up above the negative numbers (do something more clever?
#U = U/U[-1] #scale to ensure boundary condition --> makes values close to zero behave weirdly
#print(U[0], U[-1])

V_sH = U/r # equation between the equations (32) and (33)
#shift = V_Hartree[0]-V_sH[0]
#V_sH += shift #shift up (to avoid negative numbers)
#V_sH /= V_sH[0] #potential should be 1 at r=0


# PRINT TO CSV #


with open('..\\output\\A1_Task_2.csv','w') as CSV_file:
    CSV_file.write(f"Radial distance r (atomic units), theoretical Hartree potential (atomic units), calculated Hartree potential (atomic units)\n")
    for line in range(N):
            CSV_file.write(str(r[line]) + ", " + str(V_Hartree[line]) + ", " + str(V_sH[line]) + "\n")


# PLOT #
plt.plot(r, V_Hartree,          color='black', marker='', linestyle='-', label='Theoretical Hartree potential')
#plt.plot(r, electron_density,   color='black', marker='', linestyle='--', label='Theoretical electron density')
plt.plot(r, V_sH,               color='red', marker='.', linestyle='', label='Calculated Hartree potential, shifted')
#plt.plot(r, U,                  color='blue', marker='.', linestyle='', label='Solution of U = U_0 + r/r_max')
#plt.plot(r, U/r,                color='green', marker='.', linestyle='', label='U/r')

plt.xlabel('Radial distance r (atomic units)')
plt.ylabel('Energy of potential (atomic units)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
import numpy as np





def create_discretized_1D_space(start, end, number_of_points):
    a, b = start, end
    n = number_of_points

    linspace = [a+i*(b-a)/n for i in range(n+1)]
    return linspace

def Hartree_potential(radius):
    r = radius
    return (1/r) - (1 + 1/r)*np.exp(-2*r)

def create_matrix_A_task_2(number_of_points_in_discretized_1D_grid):
    return 0

print(create_discretized_1D_space(5, 10, 7))
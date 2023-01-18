import numpy as np





def create_discretized_1D_space(start, end, number_of_points):
    a, b = start, end
    n = number_of_points

    linspace = [a+i*(b-a)/n for i in range(n+1)]
    return linspace


print(create_discretized_1D_space(5, 10, 7))
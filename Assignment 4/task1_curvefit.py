import numpy as np
import matplotlib.pyplot as plt
import util

for element in ["Au", "Pt", "Rh"]:
    print(f"Fitting curve for element {element}")

    # Load potential energy as function of lattice parameter a from CSV 
    CSV = np.loadtxt(f"Assignment 4/output/TIF320_A4_T1_{element}_energy_vs_lattice_parameter_step0.01.csv", delimiter=',', skiprows=1)
    a_all = CSV[:,0]
    E_pot_all = CSV[:,1]

    # Fit second order polynomial close to the minimum
    #E_pot_min_index = np.where(E_pot_all==np.min(E_pot_all))[0][0]
    #index_to_include = 3
    a = a_all#[E_pot_min_index-index_to_include+1:E_pot_min_index+index_to_include]
    E_pot = E_pot_all#[E_pot_min_index-index_to_include+1:E_pot_min_index+index_to_include]

    fit = np.polyfit(a, E_pot, 2)
    print(f"coefficients of fit for {element}: {fit}")

    # Solution of minumum for 2nd order polynomial gives most stable lattice parameter
    a_min = -fit[1]/(2*fit[0])
    print(f"Lattice constant at minimum fitted energy for {element}: {a_min:.3f} (Å)")

    # Plot to see that is looks reasonable
    model = np.poly1d(np.polyfit(a, E_pot, 2))
    a_linspace = np.linspace(np.min(a), np.max(a), 1001, endpoint=True)
    fitted_E_pot = model(a_linspace)
    plt.plot(a, E_pot, 'x')
    plt.plot(a_linspace, fitted_E_pot, '.')
    plt.plot(a_min, np.min(fitted_E_pot), 'o')
    plt.savefig(f"Assignment 4/output/{element}_fit.pdf")
    plt.clf()
    

    # Print fitted curve to CSV
    util.print_arrays_to_CSV(f"Assignment 4/output/TIF320_A4_T1_{element}_energy_vs_lattice_parameter_step0.01_fit.csv", 
                            f"Lattice parameter (Å) for {element}", a_all, 
                            f"Potential energy (eV) for {element}", E_pot_all, 
                            f"Lattice parameter linspace (Å) for {element}", a_linspace, 
                            f"Fitted energy (eV) curve for {element}", fitted_E_pot, 
                            f"Lattice parameter (Å) for min E_pot for  {element}", [a_min],
                            f"Minimum E_pot (eV) for {element}", [np.min(fitted_E_pot)],
                            print_message=True)
import numpy as np
import matplotlib.pyplot as plt
import util

#for element in ["Au", "Pt", "Rh"]:
element = "Au"
    
# Load potential energy as function of lattice parameter a from CSV 
#CSV = np.loadtxt(f"Assignment 4/output/TIF320_A4_T1_{element}_energy_vs_lattice_parameter.csv", delimiter=',', skiprows=1)
CSV = np.loadtxt("Assignment 4/output/Au_FCC_energies.csv", delimiter=',', skiprows=1)
a = CSV[:,0]
E_pot = CSV[:,1]

# Fit second order polynomial 
fit = np.polyfit(a, E_pot, 2)
print(f"coefficients of fit: {fit}")

# Solution of minumum for 2nd order polynomial gives most stable lattice parameter
a_min = -fit[1]/(2*fit[0])
print(f"Lattice constant at minimum fitted energy: {a_min:.4f} (Å)")

# Plot to see that is looks reasonable
model = np.poly1d(np.polyfit(a, E_pot, 2))
a_linspace = np.linspace(3.5, 4.4, 1001)
fitted_E_pot = model(a_linspace)
plt.plot(a, E_pot, 'x')
plt.plot(a_linspace, fitted_E_pot, '.')
plt.plot(a_min, np.min(fitted_E_pot), 'o')
plt.savefig(f"Assignment 4/{element}_fit.pdf")

# Print fitted curve to CSV
util.print_arrays_to_CSV("Assignment 4/output/TIF320_A4_T1_Au_energy_vs_lattice_parameter_fit.csv", 
                        f"Lattice parameter for {element}", a, 
                        f"Potential energy (eV) for {element}", E_pot, 
                        f"Lattice parameter linspace for {element}", a_linspace, 
                        f"Fitted energy curve for {element}", fitted_E_pot, 
                        print_message=True)
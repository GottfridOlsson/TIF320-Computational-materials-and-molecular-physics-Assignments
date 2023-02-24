import atomic_simulation_units as asu
import numpy as np
import matplotlib.pyplot as plt
from util import print_arrays_to_CSV


# Calculated values
metals     = ["Au", "Pt", "Rh"]
E_ads_CO   = np.array([-0.124, -1.434, -1.586]) * asu.eV # Adsorption energy for CO on metals
E_ads_O    = np.array([ 0.325, -1.041, -1.751]) * asu.eV # Adsorption energy for O on metals
E_a        = -0.3 * (E_ads_O + E_ads_CO) + 0.22 * asu.eV # Activation energy for CO2 formation
A_per_site = np.array([67.994, 61.422, 57.465]) * asu.Ang**2 / 9 # Area per reaction site

# Read entropy and temperature
data_CO = np.genfromtxt("Assignment 4/output_T5/TIF320_A4_T5_entropy_vs_temperature_at_P=101325Pa_CO.csv", skip_header=1, delimiter=",")
data_O2 = np.genfromtxt("Assignment 4/output_T5/TIF320_A4_T5_entropy_vs_temperature_at_P=101325Pa_O2.csv", skip_header=1, delimiter=",")
T = data_CO[:,0] * asu.K
beta = 1 / (asu.kB * T)
assert(np.all(T == data_O2[:,0])); "Different temperatures in S vs T files!"
S_CO = data_CO[:,1] * asu.eV / asu.K
S_O2 = data_O2[:,1] * asu.eV / asu.K

# Partial pressures
p_O2 = 1 * asu.atm
p_CO = 1 * asu.atm

# Compute reaction rate for each metal
r = np.zeros((3, len(T)))
theta_CO = np.zeros((3, len(T)))
theta_O  = np.zeros((3, len(T)))
for i, metal in enumerate(metals):

    # Compute rate constants
    K_O2 = np.exp(- S_O2 / asu.kB) * np.exp( - beta * 2 * E_ads_O[i])
    K_CO = np.exp(- S_CO / asu.kB) * np.exp( - beta * 1 * E_ads_CO[i])

    # Compute fractional coverage
    theta_O[i,:] = (p_O2 * K_O2 - np.sqrt(p_O2 * K_O2) * (1 + p_CO * K_CO)) / \
            (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

    theta_CO[i,:] = (p_CO * K_CO * (np.sqrt(p_O2 * K_O2) - (1 + p_CO * K_CO))) / \
            (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

    # Compute reaction rate
    nu = 1e12 * asu.s**-1
    r_per_site = theta_O[i,:] * theta_CO[i,:] * nu * np.exp(- beta * E_a[i])

    # Scale to moles per unit area
    sites_per_area = 1 / A_per_site[i]
    print(f"{metal}: sites per area {sites_per_area / (asu.m**-2 * asu.mol)} mol / Å^2")
    r[i,:]= r_per_site * sites_per_area

# Plot
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
for i, metal in enumerate(metals):
    ax1.semilogy(T, theta_O[i,:], ["r", "g", "b"][i] + "-", label="$\\theta_{O}$, " + metal)
    ax1.semilogy(T, theta_CO[i,:], ["r", "g", "b"][i] + "--", label="$\\theta_{CO}$, " + metal)
ax1.set_ylabel("Fractional coverage")
ax1.legend()
ax1.grid()

ax2 = fig.add_subplot(2,1,2)
for i, metal in enumerate(metals):
    ax2.semilogy(T, r[i,:] / (asu.mol * asu.m**-2 * asu.s**-1), ["r", "g", "b"][i] + "-", label=metal)
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Reaction rate (mol/m$^2$s)")
ax2.grid()
ax2.legend()

plt.show()

# Export to csv
print_arrays_to_CSV(
    "Assignment 4/output_T7/rate_and_coverage_vs_temperature.csv",
    "Temperature", T / asu.K,
    "Reaction rate Au [mol / (m2 s)]", r[0,:] / (asu.mol / (asu.m**2 * asu.s)),
    "Reaction rate Pt [mol / (m2 s)]", r[1,:] / (asu.mol / (asu.m**2 * asu.s)),
    "Reaction rate Rh [mol / (m2 s)]", r[2,:] / (asu.mol / (asu.m**2 * asu.s)),
    "Fractional coverage, CO on Au", theta_CO[0,:],
    "Fractional coverage, CO on Pt", theta_CO[1,:],
    "Fractional coverage, CO on Rh", theta_CO[2,:],
    "Fractional coverage, O on Au", theta_O[0,:],
    "Fractional coverage, O on Pt", theta_O[1,:],
    "Fractional coverage, O on Rh", theta_O[2,:],
)
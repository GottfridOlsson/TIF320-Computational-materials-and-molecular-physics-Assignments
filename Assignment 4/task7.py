import atomic_simulation_units as asu
import numpy as np
import matplotlib.pyplot as plt
from util import print_arrays_to_CSV

# Calculated values
S_O2       = 198.263 * asu.J * asu.K**-1 * asu.mol**-1 # Entropy of O2 (gas)
S_CO       = 205.968 * asu.J * asu.K**-1 * asu.mol**-1 # Entropy of CO (gas)
metals     = ["Au", "Pt", "Rh"]
E_ads_CO   = np.array([-0.124, -1.434, -1.586]) * asu.eV # Adsorption energy for CO on metals
E_ads_O    = np.array([ 0.325, -1.041, -1.751]) * asu.eV # Adsorption energy for O on metals
E_a        = -0.3 * (E_ads_O + E_ads_CO) + 0.22 * asu.eV # Activation energy for CO2 formation
A_per_site = np.array([67.994, 61.422, 57.465]) * asu.Ang**2 / 9 # Area per reaction site

print(E_a / asu.eV)

# Partial pressures
p_O2 = 1 * asu.atm
p_CO = 1 * asu.atm

# Temperature range
T = np.linspace(100, 2000, 1000) * asu.K
beta = 1 / (asu.kB * T)

# Reaction rates
r = np.zeros((3, len(T)))
f = np.zeros((3, len(T)))
for i, metal in enumerate(metals):

    # Compute rate constants
    K_O2 = np.exp(- S_O2 / asu.kB) * np.exp( - beta * 2 * E_ads_O[i])
    K_CO = np.exp(- S_CO / asu.kB) * np.exp( - beta * 1 * E_ads_CO[i])

    # Compute fractional coverage
    theta_O =  (p_O2 * K_O2 - np.sqrt(p_O2 * K_O2) * (1 + p_CO * K_CO)) / \
            (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

    theta_CO = (p_CO * K_CO * (np.sqrt(p_O2 * K_O2) - (1 + p_CO * K_CO))) / \
            (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

    f[i,:] = theta_O * theta_CO

    # Compute reaction rate
    nu = 1e12 * asu.s**-1
    r_per_site = theta_O * theta_CO * nu * np.exp(- beta * E_a[i])

    # Scale to moles per unit area
    sites_per_area = 1 / A_per_site[i]
    print(f"{metal}: sites per area {sites_per_area / (asu.m**-2 * asu.mol)} mol / Å^2")
    r[i,:]= r_per_site * sites_per_area

# Plot
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.semilogy(T, np.exp(- beta * E_a[i]), "k-", label="$e^{E_a/k_bT}$")
for i, metal in enumerate(metals):
    ax1.semilogy(T, f[i,:], ["r", "g", "b"][i] + "-", label=f"$f(\\theta)$, {metal}")
ax1.set_ylabel("Fractional coverage")
ax1.legend()
ax1.grid()

ax2 = fig.add_subplot(2,1,2)
for i, metal in enumerate(metals):
    ax2.semilogy(T, r[i,:] / (asu.m**-2 * asu.s**-1), ["r", "g", "b"][i] + "-", label=metal)
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Reaction rate (mol/m$^2$s)")
ax2.grid()
ax2.legend()

plt.show()

# Export to csv
print_arrays_to_CSV(
    "Assignment 4/output_T7/r_vs_T.csv",
    "Temperature", T / asu.K,
    "Reaction rate Au [mol / (m2 s)]", r[0,:] / (asu.mol / (asu.m**2 * asu.s)),
    "Reaction rate Pt [mol / (m2 s)]", r[1,:] / (asu.mol / (asu.m**2 * asu.s)),
    "Reaction rate Rh [mol / (m2 s)]", r[2,:] / (asu.mol / (asu.m**2 * asu.s))
)
import atomic_simulation_units as asu
import numpy as np
import matplotlib.pyplot as plt

# PROPERTIES OF MOLECULES

# Masses
m_O = 15.999 * asu.u    # mass of oxygen
m_C = 12.011 * asu.u    # mass of carbon
m_O2 = m_O + m_O
m_CO = m_C + m_O

# Bond lengths
r_O2 = 1.208 * asu.Ang
r_CO = 1.128 * asu.Ang

# Moments of inertia
I_O2 = r_O2**2 * m_O / 2
I_CO = r_CO**2 * (m_O * m_C) / (m_O + m_C)

# Symmetry numbers
# "2 if the molecule has an inversion center and 1 otherwise"
sigma_O2 = 2
sigma_CO = 1

# Vibrational energy (fill in from task 5)
eps_O2 = 0.2555 * asu.eV
eps_CO = 0.2670 * asu.eV

# Energies of molecules of surface (fill in from task 4 and 6)
E_gas_O2 = 0 * asu.eV
E_gas_CO = 0 * asu.eV
E_adsorbed_O  = - 1 * asu.eV
E_adsorbed_CO = - 1 * asu.eV

# Activation energy of reaction (from PM)
E_a = - 0.3 * (E_adsorbed_O + E_adsorbed_CO) + 0.22 * asu.eV

# PROPERTIES OF SIMULATION

# Temperature range
T = np.linspace(100, 2000, 1000) * asu.K
beta = 1 / (asu.kB * T)

# Partial pressures
p_O2 = 1 * asu.atm
p_CO = 1 * asu.atm

# ENTROPIES for each degree of freedom (per particle)

# Oxygen
S_O2_trans = asu.kB * (
    np.log(
        (1 / (beta * p_O2)) * 
        (2 * np.pi * m_O2 / (beta * asu.h**2))**(3/2)
    ) 
    + 5/2
)

S_O2_rot = asu.kB * (
    np.log(8 * np.pi**2 * I_O2 / (beta * sigma_O2 * asu.h**2)) 
    + 1
)

S_O2_vib = asu.kB * (
    np.log(1 / (1 - np.exp(- beta * eps_O2))) +
    (beta * eps_O2) / (np.exp(beta * eps_O2) - 1)
)

S_O2 = S_O2_trans + S_O2_rot + S_O2_vib

# Carbon monoxide
S_CO_trans = asu.kB * (
    np.log(
        (1 / (beta * p_CO)) * 
        (2 * np.pi * m_CO / (beta * asu.h**2))**(3/2)
    ) 
    + 5/2
)

S_CO_rot = asu.kB * (
    np.log(8 * np.pi**2 * I_CO / (beta * sigma_CO * asu.h**2)) 
    + 1
)

S_CO_vib = asu.kB * (
    np.log(1 / (1 - np.exp(- beta * eps_CO))) +
    (beta * eps_CO) / (np.exp(beta * eps_CO) - 1)
)

S_CO = S_CO_trans + S_CO_rot + S_CO_vib


# Compute rate constants
K_O2 = np.exp(- S_O2 / asu.kB) * np.exp(beta * (E_gas_O2 - 2 * E_adsorbed_O))
K_CO = np.exp(- S_CO / asu.kB) * np.exp(beta * (E_gas_CO - E_adsorbed_CO))

# Compute fractional coverage
theta_O =  (p_O2 * K_O2 - np.sqrt(p_O2 * K_O2) * (1 + p_CO * K_CO)) / \
           (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

theta_CO = (p_CO * K_CO * (np.sqrt(p_O2 * K_O2) - (1 + p_CO * K_CO))) / \
           (p_O2 * K_O2 - (1 + p_CO * K_CO)**2)

# Compute reaction rate
nu = 1e12 * asu.s**-1
r = theta_O * theta_CO * nu * np.exp(- beta * E_a)

fig = plt.figure()

ax = fig.add_subplot(2,1,1)
ax.semilogy(T, theta_O, label="O")
ax.semilogy(T, theta_CO, label="CO")
ax.set_ylabel("Fractional coverage")
ax.legend()
ax.grid()

ax = fig.add_subplot(2,1,2)
ax.plot(T, r / asu.s**-1)
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Reaction rate (s^-1)")
ax.grid()


plt.show()
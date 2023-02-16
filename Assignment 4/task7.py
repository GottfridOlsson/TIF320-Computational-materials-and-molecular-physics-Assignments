import atomic_simulation_units as asu
import numpy as np

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

# Vibrational energy
eps_O2 = 0.2555 * asu.eV
eps_CO = 0.2670 * asu.eV

# Energies of molecules of surface
E_gas_O2 = 0 * asu.eV
E_gas_CO = 0 * asu.eV
E_adsorbed_O2 = - 1 * asu.eV
E_adsorbed_CO = - 1 * asu.eV

# PROPERTIES OF SIMULATION

# Temperature range
T = 25 + 273#np.linspace(100, 2000) * asu.K
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
    (beta * eps_O2) / (1 + np.exp(beta * eps_O2))
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

print(f"Entropy of O2: {S_O2 / (asu.J / (asu.K * asu.mol)):.2f} J/(K mol)")
print(f"Entropy of CO: {S_CO / (asu.J / (asu.K * asu.mol)):.2f} J/(K mol)")
print(I_O2 / (asu.kg * asu.m**2))
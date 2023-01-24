# Analytical properties of hydrogen
# source: https://en.wikipedia.org/wiki/Hydrogen_atom

import numpy as np

def ground_state_wavefunction(r, a0=1):
    return (1 / np.sqrt(np.pi)) * (1 / a0**(3/2)) * np.exp(-r / a0)

def ground_state_electron_density(r, a0=1):
    return np.exp(-2*r/a0) / (np.pi * a0**3)

def hartree_potential(r):
    return 1.0/r - (1.0 + 1.0/r) * np.exp(-2.0*r)

def ground_state_energy():
    return -0.5
# Define constants (atomic simulation units)

# Defining units
Ang     = 1                                         # Ångström
ps      = 1                                         # picosecond
eV      = 1                                         # electron volt
K       = 1                                         # kelvin

# Derived base units
m_asu   = 1 * eV * ps**2 * Ang**-2                  # asu mass unit

# SI units
s       = 1e12 * ps                                 # second
m       = 1e10 * Ang                                # meter
J       = 6.2415091e18 * eV                        # Joule
kg      = J * m**-2 * s**2                          # kilogram
N       = J * m**-1                                 # Newton
Pa      = N * m**-2                                 # Pascal

# Useful constants
u       = 1.66053907e-27 * kg                       # atomic mass unit
kB      = 1.380649e-23 * m**2 * kg * s**-2 * K**-1  # Bolzmanns constant
h       = 6.62607015e-34 * J * s                    # Plancks constant
atm     = 101325 * Pa                               # 1 atm of pressure
N_A     = 6.02214076e23                             # Avogadros constant
mol     = N_A                                       # 1 mole
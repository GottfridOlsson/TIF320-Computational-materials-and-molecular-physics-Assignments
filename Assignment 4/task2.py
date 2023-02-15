import numpy as np

DFT_energies_Au_Pt_Rh_atom = np.array([-0.158, -0.704, -1.218]) #eV/atom
T1_energies_Au_Pt_Rh_atom  = np.array([-3.146, -6.434, -7.307]) #eV/atom
cohesive_energy = DFT_energies_Au_Pt_Rh_atom - T1_energies_Au_Pt_Rh_atom

print(f"Coesive energy for Au, Pt, Rh (eV/atom): {cohesive_energy}")



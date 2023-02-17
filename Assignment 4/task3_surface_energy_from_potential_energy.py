import numpy as np

elements = ['Au', 'Pt', 'Rh']
E_bulk = [-3.146, -6.434, -7.307] # eV/atom from Task 1
E_slab = [-106.52714481388658, -195.846034, -217.5149372] # eV from Task 3
cell_lenghts = [8.860755, 8.421642, 8.145870] # angstrom
angles = np.array([60, 60, 60])*(np.pi/180) # degree
N = 27 #atoms per slab = 3*3*3

print("Area of one side of the slab")
print(f"Element     A (Å^2)")
area = []
for i, element in enumerate(elements):
    area.append(cell_lenghts[i]**2*np.sin(angles[i])) #area of parallellogram
    print(f"{element}            {area[i]:.3f}")

print("Surface energy")
print(f"Element     E_surface (eV/Å^2)")
for i, element in enumerate(elements):
    E_surface = (E_slab[i] - N*E_bulk[i]) / (2*area[i])
    print(f"{element}           {E_surface:.4f}")
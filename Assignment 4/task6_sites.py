from ase.io import write
from ase.build import add_adsorbate, molecule, fcc111

output_path_start = "Assignment 4/output_T6/"

adsorbate = molecule("CO")

# Build surfaces and add adsorbate
surface = fcc111("Rh", a=3.840, size=[3,3,3], vacuum=6.0)
add_adsorbate(surface, adsorbate, height=2, position='ontop', mol_index=-1, offset=(0,0))
add_adsorbate(surface, adsorbate, height=2, position='bridge', mol_index=-1, offset=(1,0))
add_adsorbate(surface, adsorbate, height=2, position='fcc', mol_index=-1, offset=(0,1))
add_adsorbate(surface, adsorbate, height=2, position='hcp', mol_index=-1, offset=(1,1))

# Output
write(f"{output_path_start}site_illustration_Rh.xyz", surface)
write(f"{output_path_start}site_illustration_Rh.png", surface)
write(f"{output_path_start}site_illustration_Rh_angled-view.png", surface,  rotation='10z,-80x, 5y')

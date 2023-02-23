from gpaw import GPAW, PW
from ase.io import read, write
from ase.build import fcc111
from ase.parallel import paropen


# Functions for calculations
def activation_energy(E_ads_O, E_ads_CO):
    return 0.22-0.3*(E_ads_O + E_ads_CO) #eV

def adsorption_energy_CO(E_slabWithAdsorbate_CO, E_slab, E_Adsorbate_CO):
    return E_slabWithAdsorbate_CO - E_slab - E_Adsorbate_CO

def adsorption_energy_O(E_slabWithAdsorbate_O, E_slab, E_adsorbate_O2):
    return E_slabWithAdsorbate_O - E_slab - (0.5)*E_adsorbate_O2


# Variables
surface_names = ["Au", "Pt", "Rh"]
a_from_T1 = [4.177, 3.970, 3.840] # angstrom
#E_slabs = [-78.853, -161.823, -183.081] #eV, from Task 3

adsorbate_names = ["CO", "O"]
E_adsorbate_CO  = -14.194 #eV, from Task 4
E_adsorbate_O2  = -8.728  #eV, from Task 4

positions = ['ontop', 'bridge', 'fcc', 'hcp']
output_path_start = "Assignment 4/output_T6/"

with paropen(f"{output_path_start}E_activation_CO_O.txt", 'w') as file:
    file.write(f"Surface metal, Adsorption position, Activation energy CO and O (eV)\n")

with paropen(f"{output_path_start}E_adsorption.txt", 'w') as file:
    file.write(f"Surface metal, Adsorbate, Adsorption position, Adsorption energy (eV)\n")


for i, surface_name in enumerate(surface_names):

    # 3-layered 3X3 surface cell, 6 angstrom of vacuum in +-z-direction
    surface = fcc111(surface_name, (3, 3, 3), a=a_from_T1[i], vacuum=6.0)
    calc = GPAW(xc='PBE',
            mode=PW(450),
            kpts=(4,4,1), 
            txt=f"{output_path_start}{surface_name}_GPAW_unrelaxed.txt")
    surface.set_calculator(calc)
    E_slab = surface.get_potential_energy()
    print(f"{surface_name}, E_slab={E_slab}")
    write(f"{output_path_start}{surface_name}_fcc111_initialized.xyz", surface)
    

    for position in positions:

        # Slab with adsorbate energy
        O  = read(f"{output_path_start}{surface_name}_O_{position}_relaxed.xyz")
        CO = read(f"{output_path_start}{surface_name}_CO_{position}_relaxed.xyz")
        E_slabWithAdsorbate_O  = O.get_potential_energy()
        E_slabWithAdsorbate_CO = CO.get_potential_energy()

        # Adsorption energy
        E_ads_O  = adsorption_energy_CO(E_slabWithAdsorbate_O, E_slab, E_adsorbate_O2)
        E_ads_CO = adsorption_energy_CO(E_slabWithAdsorbate_CO, E_slab, E_adsorbate_CO)

        # Activation energy
        E_act = activation_energy(E_ads_O, E_ads_CO)

        with paropen(f"{output_path_start}E_activation_CO_O.txt", 'a') as file:
            file.write(f"{surface_name}, {position}, {E_act}\n")

        with paropen(f"{output_path_start}E_adsorption.txt", 'a') as file:
            file.write(f"{surface_name}, CO, {position}, {E_ads_CO}\n")
            file.write(f"{surface_name}, O,  {position}, {E_ads_O}\n")
            
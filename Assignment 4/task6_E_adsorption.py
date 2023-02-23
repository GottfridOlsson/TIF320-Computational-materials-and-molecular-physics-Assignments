from ase.io import read

def activation_energy(E_ads_O, E_ads_CO):
    return 0.22-0.3*(E_ads_O + E_ads_CO) #eV

surface_names = ["Au", "Pt", "Rh"]
adsorbate_names = ["CO", "O"]
positions = ['ontop', 'bridge', 'fcc', 'hcp']
output_path_start = "Assignment 4/output_T6/"

with open(f"{output_path_start}E_activation_CO_O.txt", 'w') as file:
            file.write(f"Surface metal, Adsorption position, Activation energy CO and O (eV)\n")

for surface_name in surface_names:
    for position in positions:
        O  = read(f"{output_path_start}{surface_name}_O_{position}_relaxed.xyz")
        CO = read(f"{output_path_start}{surface_name}_CO_{position}_relaxed.xyz")
        
        # TODO: OBS!
        # E_ads = E_surfaceWithAdsorbant - E_surface - E_adsorbantGasForm
        # TODO: fix below
        E_ads_O  = O.get_potential_energy()
        E_ads_CO = CO.get_potential_energy()
        E_activation = activation_energy(E_ads_O, E_ads_CO)

        with open(f"{output_path_start}E_activation_CO_O.txt", 'a') as file:
            file.write(f"{surface_name}, {position}, {E_activation}\n")
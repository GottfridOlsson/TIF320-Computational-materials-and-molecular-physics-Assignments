
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

atoms = molecule('N2')
atoms.calc = EMT()
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)
potentialenergy = atoms.get_potential_energy()
vib = Vibrations(atoms)
vib.run()
vib.summary()
print(vib.summary())
#vib_energies = vib.get_energies()
quit()
thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=2, spin=0)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)
print(G)
"""


from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations

n2 = Atoms('N2', [(0, 0, 0), (0, 0, 1.1)], calculator=EMT())

BFGS(n2).run(fmax=0.01)
vib = Vibrations(n2)
vib.run()
vib.summary()
"""
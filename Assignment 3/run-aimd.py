import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.units import fs, kB
from ase.md.npt import NPT
from gpaw import GPAW

atoms = read('someStartConfiguration.xyz')
             
calc = GPAW(
    #...
    mode = 'lcao',
    xc = 'PBE',
    basis = 'dzp',
    symmetry= {'point_group ': False}, # Turn off point -group symmetry
    charge = 1, # Charged system
    txt = 'output.gpaw -out', # Redirects calculator output to this file!
)

atoms.set_calculator(calc)


dyn = NPT( # Some MD method
    #...
    atoms ,
    temperature_K = 350,
    timestep = 1000*fs, # This is not an appropriate timestep , I can tell you that!
    ttime = 20*fs, # Don’t forget the fs!
    externalstress = 0, # We don’t use the barostat , but this needs to be set anyway!
    logfile = 'mdOutput.log', # Outputs temperature (and more) to file at each timestep
    )
trajectory = Trajectory('someDynamics.traj', 'w', atoms)
dyn.attach(trajectory.write , interval =1) # Write the current positions etc. to file each timestep
dyn.run (10) # Run 10 steps of MD simulation
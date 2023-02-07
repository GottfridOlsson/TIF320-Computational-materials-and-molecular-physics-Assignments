import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.units import fs, kB
from ase.md.npt import NPT
from gpaw import GPAW

# Read snapshot where we have inserted a sodium atom into an equilibrated H2O system
atoms = read('Assignment 3/snapshots/unequilibrated-H2O-with-Na.xyz')
             
calc = GPAW(
    #...
    mode = 'lcao',
    xc = 'PBE',
    basis = 'dzp',
    symmetry= {'point_group': False}, # Turn off point -group symmetry
    charge = 1, # Charged system
    txt = 'Assignment 3/logs/output.gpaw-out', # Redirects calculator output to this file!
)

atoms.set_calculator(calc)

# NPT uses combined Nose-Hoover and Parrinello-Rahman dynamics.
dyn = NPT( # Some MD method
    #...
    atoms,
    pfactor = None,
    temperature_K = 350,
    timestep = 0.5*fs,
    ttime = 20*fs,
    externalstress = 0, # We don’t use the barostat, but this needs to be set anyway!
    logfile = 'Assignment 3/logs/mdOutput.log', # Outputs temperature (and more) to file at each timestep
    )
trajectory = Trajectory('Assignment 3/logs/nose_hoover_trajectory.traj', 'w', atoms)
dyn.attach(trajectory.write, interval=1) # Write the current positions etc. to file each timestep

# 2ps = 2000fs = 4000 * 0.5fs
timesteps = 4000
dyn.run(timesteps) # Run 10 steps of MD simulation


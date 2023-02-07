from ase.io.trajectory import Trajectory
from ase.io import write, read
traj = Trajectory('Assignment 3/Na-aimd/cluster24.traj')
atoms = traj[-1]
write('Assignment 3/snapshots/equilibrated-H2O.xyz', atoms)
write('Assignment 3/snapshots/equilibrated-H2O.png', atoms)
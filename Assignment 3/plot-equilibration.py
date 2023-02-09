import numpy as np
from util import *
import matplotlib.pyplot as plt

# Read data from log file
data = np.genfromtxt('Assignment 3/Na-aimd/cluster24.log', skip_header=1)
time = data[:,0]
Etot = data[:,1]
Epot = data[:,2]
Ekin = data[:,3]
T    = data[:,4]

T_instant = T
T_avg     = np.cumsum(T) / np.cumsum(np.ones_like(T))
T_target  = 350 * np.ones_like(time)

#plt.plot(time, T_instant,   label="Instantaneous temperature")
#plt.plot(time, T_avg,       label="Running average")
#plt.plot(time, T_target,    "k--", label=f"T = {T_target[0]}K")
plt.plot(time, Etot,         "k", label=f"E_tot")
plt.grid()
plt.legend()
plt.xlim(0, 3)
plt.xlabel("Time (ps)")
plt.ylabel("Temperature (K)")
plt.savefig('Assignment 3/Na-aimd/cluster24-simulation-temperature-vs-time.pdf')

print_arrays_to_CSV('Assignment 3/Na-aimd/their-simulation-without-Na-temperature-energy-vs-time.csv',
    'Time (ps)', time,
    'Instantaneous temperature (K)', T_instant,
    'Running average temperature (K)', T_avg,
    'Target temperature (K)', T_target, 
    'Total energy (eV)', Etot, 
    'Potential energy (eV)', Epot, 
    'Kinetic energy (eV)', Ekin
)
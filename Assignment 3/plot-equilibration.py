import numpy as np
from util import *
import matplotlib.pyplot as plt

# Read data from log file
data = np.genfromtxt('Assignment 3/logs/mdOutput.log', skip_header=1)
time = data[:,0]
Etot = data[:,1]
Epot = data[:,2]
Ekin = data[:,3]
T    = data[:,4]

T_instant = T
T_avg     = np.cumsum(T) / np.cumsum(np.ones_like(T))
T_target  = 345 * np.ones_like(time)

plt.plot(time, T_instant,   label="Instantaneous temperature")
plt.plot(time, T_avg,       label="Running average")
plt.plot(time, T_target,    "k--", label="T = 350K")
plt.grid()
plt.legend()
plt.ylim(200, 500)
plt.xlabel("Time (ps)")
plt.ylabel("Temperature (K)")
plt.savefig('Assignment 3/logs/our-simulation-temperature-vs-time.pdf')

print_arrays_to_CSV('Assignment 3/logs/our-simulation-temperature-vs-time.csv',
    'Time (ps)', time,
    'Instantaneous temperature (K)', T_instant,
    'Running average temperature (K)', T_avg,
    'Target temperature (K)', T_target
)
#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

trajectories = np.genfromtxt('./out/trajectories.csv', delimiter=',', skip_header=1)

plt.figure(figsize=(10,4))
plt.plot(trajectories[:, 0], trajectories[:, 1])
#plt.xlim(0.0,0.1e9)
#plt.ylim(0.0,3000)
plt.savefig('./plots/positions.pdf')

plt.figure(figsize=(10,4))
plt.plot(trajectories[:, 0], trajectories[:, 2])
#plt.xlim(2200,2300)
plt.savefig('./plots/velocities.pdf')

plt.figure(figsize=(10,4))
plt.plot(trajectories[:, 0], trajectories[:, 3])
#plt.xlim(2200,2300)
#plt.ylim(-3.0e-11,3.0e-11)
plt.savefig('./plots/accelerations.pdf')

#local_energy = np.genfromtxt('./out/local_energy.csv', delimiter=',', skip_header=1)
#fig, ax = plt.subplots(figsize=(10,4))
#ax.plot(local_energy[:, 0], local_energy[:, 1])
#fig.savefig('./plots/local_energy.pdf')

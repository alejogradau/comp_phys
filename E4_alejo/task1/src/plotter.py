#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

font = {'size'   : 16}
mpl.rc('font', **font)

equilibrated = True

trajectories = np.genfromtxt('./out/trajectories.csv', delimiter=',', skip_header=1)
equilibration = 10  # 10 ms

if(equilibrated):
    plt.figure(figsize=(8,8))
    plt.scatter(trajectories[:, 0], trajectories[:, 1], s=0.5)
    plt.plot(trajectories[:, 0], trajectories[:, 1])
    plt.xlim(0.0,2.0)
    plt.xlabel("Time (ms)")
    plt.ylabel("Position (nm)")
    plt.ylim(-80,80)
    plt.savefig('./plots/positions.pdf')

    plt.figure(figsize=(8,8))
    plt.scatter(trajectories[:, 0], trajectories[:, 2], s=0.5)
    plt.plot(trajectories[:, 0], trajectories[:, 2])
    plt.xlim(0.0,2.0)
    plt.xlabel("Time (ms)")
    plt.ylabel("Velocity (mm/s)")
    plt.ylim(-2,2)
    plt.savefig('./plots/velocities.pdf')

    #plt.figure(figsize=(10,4))
    #plt.plot(trajectories[:, 0], trajectories[:, 3])
    #plt.xlim(2200,2300)
    #plt.ylim(-3.0e-11,3.0e-11)
    #plt.savefig('./plots/accelerations.pdf')

else:
    plt.figure(figsize=(8,8))
    plt.scatter(trajectories[:, 0] - equilibration, trajectories[:, 1], s=0.5)
    plt.plot(trajectories[:, 0] - equilibration, trajectories[:, 1])
    plt.xlim(0.0,2.0)
    plt.xlabel("Time (ms)")
    plt.ylabel("Position (nm)")
    plt.ylim(-80,80)
    plt.savefig('./plots/positions.pdf')

    plt.figure(figsize=(8,8))
    plt.scatter(trajectories[:, 0] - equilibration, trajectories[:, 2], s=0.5)
    plt.plot(trajectories[:, 0] - equilibration, trajectories[:, 2])
    plt.xlim(0.0,2.0)
    plt.xlabel("Time (ms)")
    plt.ylabel("Velocity (mm/s)")
    plt.ylim(-2,2)
    plt.savefig('./plots/velocities.pdf')

    #plt.figure(figsize=(10,4))
    #plt.plot(trajectories[:, 0], trajectories[:, 3])
    #plt.xlim(2200,2300)
    #plt.ylim(-3.0e-11,3.0e-11)
    #plt.savefig('./plots/accelerations.pdf')

#local_energy = np.genfromtxt('./out/local_energy.csv', delimiter=',', skip_header=1)
#fig, ax = plt.subplots(figsize=(10,4))
#ax.plot(local_energy[:, 0], local_energy[:, 1])
#fig.savefig('./plots/local_energy.pdf')

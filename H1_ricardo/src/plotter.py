#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

energy = True
temperature = True


if energy:
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
#    ax.set_title('FCC Lattice')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    plt.legend()

    fig.savefig('./output/kinetic.pdf')

if temperature:
    array = np.genfromtxt('./output/temperature.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('FCC Lattice')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend()

    fig.savefig('./output/temperature.pdf')

#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

fig1_energy_volume = True

if fig1_energy_volume:
    array = np.genfromtxt('csv/energy_volume.csv',
    delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(8,4))
    ax.scatter(array[:, 0], array[:, 1], s=4)

    ax.set_xlabel('Volume (Å^3)')
    ax.set_ylabel('Energy (eV/unit cell)')

    fig.savefig('plots/fig1_energy_volume.pdf')

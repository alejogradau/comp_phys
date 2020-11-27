#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

energy = True


if energy:
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic E")
    ax.plot(array[:, 0], array[:, 2], label="Potential E")
    ax.plot(array[:, 0], array[:, 3], label="Total E")
    ax.set_title('N-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    plt.legend()

    fig.savefig('./output/kinetic.pdf')

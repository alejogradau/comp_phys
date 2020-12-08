#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

local_energy = np.genfromtxt('./out/local_energy.csv', delimiter=',', skip_header=1)
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(local_energy[:, 0], local_energy[:, 1])
fig.savefig('./plots/local_energy.pdf')

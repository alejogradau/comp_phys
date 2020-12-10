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

signal_50 = np.genfromtxt('./out/powerspectrum_50dt.csv', delimiter=',', skip_header=1)
signal_25 = np.genfromtxt('./out/powerspectrum_25dt.csv', delimiter=',', skip_header=1)


plt.figure(figsize=(8,8))
plt.plot(signal_50[:, 0], signal_50[:, 1], label=r'$d\tau = 25*dt$')
plt.plot(signal_25[:, 0], signal_25[:, 1], label=r'$d\tau = 50*dt$')
plt.xlim(0.0,10.0)
plt.xlabel("Frequency (kHz)")
plt.ylabel("Power (arb. units)")
#plt.ylim(-80,80)
plt.legend()
plt.savefig('./plots/powerspectrum.pdf')

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

#signal = np.genfromtxt('./out/powerspectrum_dt.csv', delimiter=',', skip_header=1)
#signal_50 = np.genfromtxt('./out/powerspectrum_50dt.csv', delimiter=',', skip_header=1)
#signal_25 = np.genfromtxt('./out/powerspectrum_25dt.csv', delimiter=',', skip_header=1)


factor = 25

signal_low = np.genfromtxt('./out/powerspectrum_low_'
            + str(factor) + 'dt.csv', delimiter=',', skip_header=1)
signal_high = np.genfromtxt('./out/powerspectrum_high_'
            + str(factor) + 'dt.csv', delimiter=',', skip_header=1)


plt.figure(figsize=(8,8))
plt.plot(signal_low[:, 0], signal_low[:, 1], label=r'Case Low: $d\tau ='
         + str(factor) + ' dt$')
plt.plot(signal_high[:, 0], signal_high[:, 1], label=r'Case High: $d\tau ='
         + str(factor) + ' dt$')
#plt.plot(signal[:, 0], signal[:, 1], label=r'$d\tau = dt$')
#plt.plot(signal_25[:, 0], signal_25[:, 1], label=r'$d\tau = 25*dt$')
#plt.plot(signal_50[:, 0], signal_50[:, 1], label=r'$d\tau = 50*dt$')
plt.xlim(0.0,10.0)
#plt.ylim(-80,80)
plt.xlabel("Frequency (kHz)")
plt.ylabel("Power (arb. units)")
plt.legend()
plt.savefig('./plots/powerspectrum_' + str(factor) + 'dt.pdf')

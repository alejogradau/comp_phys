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

zoom = 1000  # microseconds
factor = 50

autocorrelation_low = np.genfromtxt('./out/autocorrelation_low.csv',
                                    delimiter=',', skip_header=1)
autocorrelation_high = np.genfromtxt('./out/autocorrelation_high.csv',
                                    delimiter=',', skip_header=1)

#len = len(autocorrelation_low)
#len = m.floor(len/factor)
#power_low = np.zeros(len)
#power_high = np.zeros(len)

signal_low = np.genfromtxt('./out/powerspectrum_low_'
            + str(factor) + 'dt.csv', delimiter=',', skip_header=1)
signal_high = np.genfromtxt('./out/powerspectrum_high_'
            + str(factor) + 'dt.csv', delimiter=',', skip_header=1)

power_low = signal_low[:, 1]
power_high = signal_high[:, 1]
frequencies_high = signal_high[:, 0]
frequencies_low = signal_low[:, 0]

#  PLOT AUTOCORRELATION
plt.figure(figsize=(10,5))
plt.plot(autocorrelation_low[0:zoom, 0], autocorrelation_low[0:zoom, 1], label='Low')
plt.plot(autocorrelation_high[0:zoom, 0], autocorrelation_high[0:zoom, 1], label='High')
plt.xlabel("Time (ms)")
plt.ylabel("Normalized VACF")
plt.legend()
plt.savefig('./plots/autocorrelation.pdf')

#PLOT FFT
plt.figure(figsize=(8,8))
plt.plot(frequencies_low, power_low, label=r' Case Low: $d\tau = '
                                        + str(factor) + '*dt$')
plt.plot(frequencies_high, power_high, label=r' Case High: $d\tau = '
                                        + str(factor) + '*dt$')
plt.xlim(0.0,10.0)
plt.xlabel("Frequency (kHz)")
plt.ylabel("Power (arb. units)")
#plt.ylim(-80,80)
plt.legend()
plt.savefig('./plots/powerspectrum_' + str(factor) + 'dt.pdf')

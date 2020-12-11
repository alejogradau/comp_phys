#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math as m

font = {'size'   : 16}
mpl.rc('font', **font)

M = 100
factor = 50
len = 99900
len = m.floor(len/factor)
power_low = np.zeros(len)
power_high = np.zeros(len)

#len_1 = 99900
#len_25 = m.floor(len_1/25)
#len_50 = m.floor(len_1/50)

#power_1 = np.zeros(len_1)
#power_25 = np.zeros(len_25)
#power_50 = np.zeros(len_50)

#for i in range(M):
#    signal_1 = np.genfromtxt('./out/powerspectrum_'
#                + str(i) + '_1dt.csv', delimiter=',', skip_header=1)
#    power_1 = power_1 + signal_1[:, 1]
#
#power_1 /= M
#frequencies_1 = signal_1[:, 0]

for i in range(M):
    signal_low = np.genfromtxt('./out/powerspectrum_low_' + str(i) + '_'
                + str(factor) + 'dt.csv', delimiter=',', skip_header=1)
    signal_high = np.genfromtxt('./out/powerspectrum_high_' + str(i) + '_'
                + str(factor) + 'dt.csv', delimiter=',', skip_header=1)
    power_low = power_low + signal_low[:, 1]
    power_high = power_high + signal_high[:, 1]

power_low /= M
power_high /= M
frequencies = signal_low[:, 0]

plt.figure(figsize=(8,8))
plt.plot(frequencies, power_low, label=r' Case Low: $d\tau = '
                                        + str(factor) + '*dt$')
plt.plot(frequencies, power_high, label=r' Case High: $d\tau = '
                                        + str(factor) + '*dt$')
plt.xlim(0.0,10.0)
plt.xlabel("Frequency (kHz)")
plt.ylabel("Power (arb. units)")
#plt.ylim(-80,80)
plt.legend()
plt.savefig('./plots/powerspectrum_50dt_M' + str(M) + '.pdf')

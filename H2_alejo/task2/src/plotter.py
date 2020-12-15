#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

Z = 'Z = 1.7'
Zname = '_Z1.7'
fname1 = 'radial_density'
fname2 = 'histogram'
figsname = 'histogram_density_comparison'
datapath = './out/'
figpath = './plots/'
dataext = '.csv'
figext = '.pdf'

font = {'size'   : 16}
mpl.rc('font', **font)

radial_density = np.genfromtxt(datapath + fname1 + Zname + dataext,
                                delimiter=',', skip_header=1)
histogram = np.genfromtxt(datapath + fname2 + dataext,
                            delimiter=',', skip_header=1)

x = histogram[:, 0]
y1 = histogram[:, 1]
y2 = radial_density[:, 1]

#plt.figure(figsize=(10,4))
#plt.scatter(radial_density[:, 0]+1, radial_density[:, 1], s=5)
#plt.plot(radial_density[:, 0]+1, radial_density[:, 1], label = 'Z=2')
#plt.xlim(2200,2300)
#plt.savefig(figpath + fname1 + figext)

#plt.figure(figsize=(10,4))
#plt.scatter(histogram[:, 0], histogram[:, 1], s=5)
#plt.plot(histogram[:, 0], histogram[:, 1])
#plt.xlim(2200,2300)
#plt.savefig(figpath + fname2 + figext)

fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(111)
ax1.plot(x, y1, 'b-')
ax1.set_ylabel('Histogram', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(x, y2, 'r-')
ax2.set_ylabel(r'Radial probability density $\rho(r)$', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

fig.suptitle(Z, fontsize=16)
plt.savefig(figpath + figsname + Zname + figext)

#local_energy = np.genfromtxt('./out/local_energy.csv', delimiter=',', skip_header=1)
#fig, ax = plt.subplots(figsize=(10,4))
#ax.plot(local_energy[:, 0], local_energy[:, 1])
#fig.savefig('./plots/local_energy.pdf')

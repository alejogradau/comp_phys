#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

datapath = './out/task2/'
figpath = './plots/'
dataext = '.csv'
figext = '.pdf'

Z = 2.0
title = 'Z = ' + str(Z)
Zname = '_Z' + str(Z)
fname1 = 'radial_density'
fname2 = 'histogram'
figsname1 = 'histogram_density_comparison'

fname3 = 'local_energy_alpha0.10_run0'
figsname2 = fname3

font = {'size'   : 16}
mpl.rc('font', **font)

histogram_density_comparison = False
local_energy = True

if histogram_density_comparison:
    radial_density = np.genfromtxt(datapath + fname1 + Zname + dataext,
                                delimiter=',', skip_header=1)
    histogram = np.genfromtxt(datapath + fname2 + dataext,
                            delimiter=',', skip_header=1)

    x_radial = np.zeros(21)
    dr = 0.2364817
    for i in range(21):
        x_radial[i] = dr*(i+0.5)
    x = x_radial
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
    ax1.set_xlabel('Radial distance (Å)')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(x, y2, 'r-')
    ax2.set_ylabel(r'Radial probability density $\rho(r)$', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    fig.suptitle(title, fontsize=16)
    plt.savefig(figpath + figsname1 + Zname + figext)

if local_energy:
    local_energy = np.genfromtxt(datapath + fname3 + dataext,
                                 delimiter=',', skip_header=1)
    fig, ax = plt.subplots(figsize=(10,6))
    ax.scatter(local_energy[:, 0], local_energy[:, 1])
    ax.set_xscale('log')
    #plt.xlim(1,100)
    ax.set_xlabel('Configuration number')
    ax.set_ylabel('Local energy (eV)')
    fig.savefig(figpath + figsname2 + figext)

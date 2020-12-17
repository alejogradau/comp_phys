#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

datapath = './out/'
figpath = './plots/'
dataext = '.csv'
figext = '.pdf'

Z = 1.7
title = 'Z = ' + str(Z)
Zname = '_Z' + str(Z)
fname1 = 'radial_density'
fname2 = 'histogram_Z' + str(Z)
figsname1 = 'histogram_density_comparison'

fname3 = 'local_energy_alpha0.10_run'
figsname2 = fname3

font = {'size'   : 16}
mpl.rc('font', **font)

histogram_density_comparison = False
local_energy = True
autoc = True
block = True

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
    local_energy1 = np.genfromtxt(datapath + fname3 + "0" + dataext, delimiter=',', skip_header=1)
#    local_energy2 = np.genfromtxt(datapath + fname3 + "1" + dataext, delimiter=',', skip_header=1)
#    local_energy3 = np.genfromtxt(datapath + fname3 + "2" + dataext, delimiter=',', skip_header=1)
#    local_energy4 = np.genfromtxt(datapath + fname3 + "3" + dataext, delimiter=',', skip_header=1)
#    local_energy5 = np.genfromtxt(datapath + fname3 + "4" + dataext, delimiter=',', skip_header=1)
                    
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(local_energy1)
    ax.set_xscale('log')
#    plt.xlim(1,10000)
    ax.set_xlabel('Configuration number')
    ax.set_ylabel('Local energy (eV)')
    fig.savefig(datapath + figsname2 + figext)

if autoc:
    autocorrelation = np.genfromtxt('./out/autocorrelation_local_energy_alpha0.10_run0.csv', delimiter=',', skip_header=1)
    fig, ax = plt.subplots(figsize=(10,6))
    ax.set_xscale('log')
    ax.plot(autocorrelation[:, 0], autocorrelation[:, 1])
    ax.set_xlabel('k')
    ax.set_ylabel('ACF')
    fig.savefig('./out/autocorrelation.pdf')
    
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(autocorrelation[0:400, 0], autocorrelation[0:400, 1])
    ax.set_xlabel('k')
    ax.set_ylabel('ACF')
    fig.savefig('./out/autocorrelation_400.pdf')


if block:
    block_averaging = np.genfromtxt('./out/block_averaging_local_energy_alpha0.10_run0.csv', delimiter=',', skip_header=1)
    fig, ax = plt.subplots(figsize=(10,6))
    ax.set_xscale('log')
    ax.set_xlabel('Block Size')
    ax.set_ylabel('Statistical Inefficiency')
    plt.scatter(block_averaging[:, 0], block_averaging[:, 1], s=5)
#    ax.plot(autocorrelation[0:400, 0], autocorrelation[0:400, 1])
    fig.savefig('./out/block_averaging.pdf')

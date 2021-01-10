#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import glob

datapath = './out/'
figpath = './plots/'
dataext = '.csv'
figext = '.pdf'

fname_grid = 'grid'
fname_prob = 'prob_density'
figname_prob = 'prob_density'

suffix_x = '_x'
suffix_p = '_p'

font = {'size'   : 16}
mpl.rc('font', **font)

prob_density_x = True
prob_density_p = True

if prob_density_x:

    data_x = np.genfromtxt(datapath + fname_grid + suffix_x + dataext,
                                delimiter=',', skip_header=1)

    data_y = np.genfromtxt(datapath + fname_prob + suffix_x + dataext,
                                delimiter=',', skip_header=1)

    grid = data_x[:,1]
    prob = data_y[:,1]

    plt.figure(figsize=(10,6))
    plt.plot(grid, prob, '--', c='k', label='Gaussian wave packet \n position space')
    plt.xlabel(r'x (Å)')
    plt.ylabel(r'Probability density distribution')
    plt.legend()
    plt.savefig(figpath + figname_prob + suffix_x + figext)

if prob_density_p:

    data_p = np.genfromtxt(datapath + fname_grid + suffix_p + dataext,
                                delimiter=',', skip_header=1)

    data_y = np.genfromtxt(datapath + fname_prob + suffix_p + dataext,
                                delimiter=',', skip_header=1)

    grid = data_p[:,1]
    prob = data_y[:,1]

    plt.figure(figsize=(10,6))
    plt.plot(grid, prob, '--', c='k', label='Gaussian wave packet \n momentum space')
    plt.xlabel(r'Momentum (eV * fs / Å )')
    #plt.xlim(-5,5)
    plt.ylabel(r'Probability density distribution')
    plt.legend()
    plt.savefig(figpath + figname_prob + suffix_p + figext)

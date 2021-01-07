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
    plt.xlabel(r'x (Å)')
    plt.ylabel(r'Probability density distribution')
    plt.legend()
    plt.savefig(figpath + figname_prob + suffix_p + figext)

if other:

    alpha = np.zeros(n_alphas)
    mean_energy = np.zeros(n_alphas)
    counter = 0
    for i in range(n_alphas):
        alpha[i] = 0.06 + i*dalpha
    for fname in filelist:
        alpha_energy = np.genfromtxt(fname, delimiter=',', skip_header=1)
        mean_energy[counter] = alpha_energy[0]
        print(alpha[counter], mean_energy[counter])
        counter +=1

    plt.figure(figsize=(10,6))
    plt.plot(alpha, mean_energy)
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'Mean energy (eV)')
    plt.savefig('alpha_vs_energy.pdf')

if steepest:

    fig, ax = plt.subplots(figsize=(10,6))
    for fname in filelist:
        alpha_trajectory = np.genfromtxt(fname, delimiter=',', skip_header=1)
        p = alpha_trajectory[:, 0]
        alpha = alpha_trajectory[:, 1]
        ax.plot(p, alpha, label=fname)
        #plt.xlim(1,100)
    plt.legend()
    ax.set_xlabel(r'Steepest descent step $p$')
    ax.set_ylabel(r'$\alpha$ ($Å^{-1}$)')
    fig.savefig('steepest_descent.pdf')


if nablas:
    nablas = np.genfromtxt(datapath + fname4 + dataext,
                                delimiter=',', skip_header=1)
    x = nablas[:, 0]
    nabla_i = nablas[:, 1]
    nabla_mean = nablas[:, 2]

    fig = plt.figure(figsize=(20,6))
    ax1 = fig.add_subplot(111)
    ax1.plot(x, nabla_i)
    ax1.set_ylabel('nabla_i', color='b')
    ax1.set_xlabel('Configuration number')
    #ax1.set_xlim(284530,284570)
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(x, nabla_mean, 'r-')
    ax2.set_ylabel('nabla_mean', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    fig.suptitle(fname4, fontsize=16)
    plt.savefig(figpath + fname4 + figext)

if histogram_density_comparison:
    radial_density = np.genfromtxt(datapath + fname1 + Zname + dataext,
                                delimiter=',', skip_header=1)
    histogram = np.genfromtxt(datapath + fname2 + dataext,
                            delimiter=',', skip_header=1)

    x_radial = np.zeros(21)
    density = Z_cb * 4 * rad_sq * exp(-2 * Z * rad);
    dr = 0.242800  # bin_size: to be changed for each run!
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

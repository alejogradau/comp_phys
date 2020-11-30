#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

# skip_header skips the first
# row in data.csv

energy = False
energy_shifted = True
temperature = True
pressure = True
a0_ev = False

font = {'size'   : 16}
mpl.rc('font', **font)


if energy:
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    fig = plt.figure(figsize=(10,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

    ax0 = plt.subplot(gs[0])
    line0, = ax0.plot(array[:, 0], array[:, 1], label="Kinetic Energy", color="C0")

    ax1 = plt.subplot(gs[1], sharex = ax0)
    line1, = ax1.plot(array[:, 0], array[:, 2], label="Potential Energy", color="C1")
    line2, = ax1.plot(array[:, 0], array[:, 3], label="Total Energy", color="C2")

    plt.setp(ax0.get_xticklabels(), visible=False)
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)

    ax0.legend((line0, line1, line2), ("Kinetic Energy", "Potential Energy", "Total Energy"), loc="lower right", fontsize='small')
    ax1.set_xlabel('Time (ps)')
    ax0.set_ylabel('Energy (eV)')
    plt.subplots_adjust(hspace=.05)

    fig.savefig('./plots/energy.pdf')

if energy_shifted:
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
#    ax.set_title('N-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    ax.set_ylim(-50,50)
    plt.legend(loc="lower right", fontsize='small')

    fig.savefig('./plots/energy.pdf')


if temperature:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/temperature.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="Instantaneous temperature")
    #ax.set_title('Time evolution of the instantaneous temperature')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (K)')
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/temperature.pdf')

if pressure:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/pressure.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="Instantaneous pressure")
    #ax.set_title('Time evolution of the instantaneous pressure')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    ax.set_ylim(-50,50)
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/pressure.pdf')

if a0_ev:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/a0.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1])
    #ax.set_title('Time evolution of the instantaneous pressure')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Lattice parameter (Å)')
    #ax.set_ylim(-50,50)
    #plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/a0.pdf')

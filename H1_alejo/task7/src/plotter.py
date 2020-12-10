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

energy = True
energy_scaled = False
temperature = True
pressure = True
lattice_param = True
positions = True

font = {'size'   : 16}
mpl.rc('font', **font)

if energy:
    array = np.genfromtxt('./out/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
#    ax.set_title('N-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    plt.legend(loc="lower right", fontsize='small')

    fig.savefig('./plots/energy.pdf')

if energy_scaled:
    array = np.genfromtxt('./out/energy.csv', delimiter=',', skip_header=1)

    fig = plt.figure()
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
    plt.subplots_adjust(hspace=.0)

    fig.savefig('./plots/energy_scaled.pdf')

if temperature:
    array = np.genfromtxt('./out/temperature.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Instantaneous temperature")
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (K)')
    plt.legend(loc="lower right", fontsize='small')

    fig.savefig('./plots/temperature.pdf')

if pressure:
    array = np.genfromtxt('./out/pressure.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Instantaneous pressure")
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend(loc="lower right", fontsize='small')

    fig.savefig('./plots/pressure.pdf')

if lattice_param:
    array = np.genfromtxt('./out/a0.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1])
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Lattice Parameter (Å)')
    #plt.legend()

    fig.savefig('./plots/a0.pdf')


if positions:
    array = np.genfromtxt('./out/positions.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label = 'x,y,z-coord. particle 0', color="C0")
    ax.plot(array[:, 0], array[:, 2], color="C0")
    ax.plot(array[:, 0], array[:, 3], color="C0")
    ax.plot(array[:, 0], array[:, 4], label = 'x,y,z-coord. particle 127', color="C1")
    ax.plot(array[:, 0], array[:, 5], color="C1")
    ax.plot(array[:, 0], array[:, 6], color="C1")
    ax.plot(array[:, 0], array[:, 7], label = 'x,y,z-coord. particle 255', color="C2")
    ax.plot(array[:, 0], array[:, 8], color="C2")
    ax.plot(array[:, 0], array[:, 9], color="C2")
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Position (Å)')
    plt.legend(loc="lower right", fontsize='small')

    fig.savefig('./plots/positions.pdf')
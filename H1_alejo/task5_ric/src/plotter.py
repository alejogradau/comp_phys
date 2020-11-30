#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# skip_header skips the first
# row in data.csv

energy = True
energy_scaled = True
temperature = True
pressure = True
temperature_avg = True
pressure_avg = True
lattice_param = True
positions = True

if energy:
    array = np.genfromtxt('./out/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
#    ax.set_title('N-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    plt.legend()

    fig.savefig('./out/energy.pdf')
    
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
    
    ax0.legend((line0, line1, line2), ("Kinetic Energy", "Potential Energy", "Total Energy"), loc="upper right")
    ax1.set_xlabel('Time (ps)')
    ax0.set_ylabel('Energy (eV)')
    plt.subplots_adjust(hspace=.0)

    fig.savefig('./out/energy_scaled.pdf')

if temperature:
    array = np.genfromtxt('./out/temperature.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend()

    fig.savefig('./out/temperature.pdf')

if pressure:
    array = np.genfromtxt('./out/pressure.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend()

    fig.savefig('./out/pressure.pdf')

if temperature_avg:
    array = np.genfromtxt('./out/temperature_avg.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend()

    fig.savefig('./out/temperature_avg.pdf')

if pressure_avg:
    array = np.genfromtxt('./out/pressure_avg.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend()

    fig.savefig('./out/pressure_avg.pdf')

if lattice_param:
    array = np.genfromtxt('./out/a0.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
#    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Lattice Parameter (Å)')
    plt.legend()

    fig.savefig('./out/a0.pdf')


if positions:
    array = np.genfromtxt('./out/positions.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1])
    ax.plot(array[:, 0], array[:, 2])
    ax.plot(array[:, 0], array[:, 3])
    ax.plot(array[:, 0], array[:, 4])
    ax.plot(array[:, 0], array[:, 5])
    ax.plot(array[:, 0], array[:, 6])
    ax.plot(array[:, 0], array[:, 7])
    ax.plot(array[:, 0], array[:, 8])
    ax.plot(array[:, 0], array[:, 9])
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Position (Å)')
    plt.legend()

    fig.savefig('./out/positions.pdf')

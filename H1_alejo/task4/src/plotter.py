#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

energy = True
temperature = True
pressure = True
temperature_avg = False
pressure_avg = False

font = {'size'   : 16}
mpl.rc('font', **font)

if energy:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
    ax.set_title('Energy evolution during equilibration at T=700K')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    ax.set_ylim(-35,35)
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/energy.pdf')

if temperature:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/temperature.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="Instantaneous temperature")
    ax.set_title('Time evolution of the instantaneous temperature')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (K)')
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/temperature.pdf')

if pressure:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/pressure.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="Instantaneous pressure")
    ax.set_title('Time evolution of the instantaneous pressure')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    ax.set_ylim(-50,50)
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/pressure.pdf')

if temperature_avg:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/temperature_avg.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/temperature_avg.pdf')

if pressure_avg:
    fig, ax = plt.subplots(figsize=(10,6))
    array = np.genfromtxt('./output/pressure_avg.csv', delimiter=',', skip_header=1)

    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend(fontsize='small', loc='lower right')

    fig.savefig('./plots/pressure_avg.pdf')

#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

energy = True
temperature = True
pressure = True
temperature_avg = True
pressure_avg = True

if energy:
    array = np.genfromtxt('./output/energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="Kinetic Energy")
    ax.plot(array[:, 0], array[:, 2], label="Potential Energy")
    ax.plot(array[:, 0], array[:, 3], label="Total Energy")
    ax.set_title('N-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    plt.legend()

    fig.savefig('./output/energy.pdf')

if temperature:
    array = np.genfromtxt('./output/temperature.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend()

    fig.savefig('./output/temperature.pdf')

if pressure:
    array = np.genfromtxt('./output/pressure.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend()

    fig.savefig('./output/pressure.pdf')

if temperature_avg:
    array = np.genfromtxt('./output/temperature_avg.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Temperature (ºK)')
    plt.legend()

    fig.savefig('./output/temperature_avg.pdf')

if pressure_avg:
    array = np.genfromtxt('./output/pressure_avg.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label="")
    ax.set_title('N-particle system Time Average')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Pressure (bar)')
    plt.legend()

    fig.savefig('./output/pressure_avg.pdf')

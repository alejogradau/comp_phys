#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv

spectrum_code1 = False
spectrum_code4 = True
spectrums_sum_code4 = False
signal = False
positions = False
abs_displacements = False
energy = False
kinetic = False
kinetic_spectrum = False


if spectrum_code1:
    array = np.genfromtxt('powerspectrum_code1.csv',
    delimiter=',', skip_header=1)

    fig, ax = plt.subplots()
    ax.plot(array[:, 0], array[:, 1])

    ax.set_xlabel('frequency (arb.unit)')
    ax.set_ylabel('Power spectrum (arb.unit)')
    ax.grid()

    fig.savefig('powerspectrum_code1.pdf')

if spectrum_code4:
    array = np.genfromtxt('powerspectrum_code4.csv',
    delimiter=',', skip_header=1)

    fig, ax = plt.subplots()
    ax.set_xlim(-100,100)
    ax.plot(array[:, 0], array[:, 1])

    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Power Spectrum (Å^-2)')
    ax.grid()

    fig.savefig('powerspectrum_code4.pdf')

if spectrums_sum_code4:
    q1 = np.genfromtxt('powerspectrum_q1.csv',
    delimiter=',', skip_header=1)
    q2 = np.genfromtxt('powerspectrum_q2.csv',
    delimiter=',', skip_header=1)
    q3 = np.genfromtxt('powerspectrum_q3.csv',
    delimiter=',', skip_header=1)

    frequencies = q1[:, 0]

    spec_q1 = q1[:, 1]
    spec_q2 = q2[:, 1]
    spec_q3 = q3[:, 1]
    spec_sum = spec_q1 + spec_q2 + spec_q3

    fig, ax = plt.subplots()
    ax.set_xlim(-100,100)
    ax.plot(frequencies, spec_sum)

    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Power Spectrum (Å^-2)')
    ax.grid()

    fig.savefig('powerspectrum_sum.pdf')

if signal:
    array = np.genfromtxt('signal.csv', delimiter=',',
    skip_header=1)

    fig, ax = plt.subplots()
    ax.plot(array[:, 0], array[:, 1])

    ax.set_xlabel('time (arb.unit)')
    ax.set_ylabel('signal (arb.unit)')
    ax.grid()

    fig.savefig('signal.pdf')


if positions:
    array = np.genfromtxt('positions.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1],
            array[:, 0], array[:, 2],
            array[:, 0], array[:, 3])
    ax.set_title('Three-particle system time-dependent displacements')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Displacement (Å)')
    ax.set_xlim(0.0,0.25)

    fig.savefig('displacement.pdf')

if abs_displacements:
    array = np.genfromtxt('abs_total_displacements.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1])
    ax.set_title('Three-particle system time-dependent absolute displacements')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Displacement (Å)')
    ax.set_xlim(0.0,0.25)

    fig.savefig('absolute_displacements.pdf')


if energy:
    array = np.genfromtxt('energy.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1], label='Kinetic')
    ax.plot(array[:, 0], array[:, 2], label='Potential')
    ax.plot(array[:, 0], array[:, 3], label='Total Energy')
    ax.set_title('Three-particle system time-dependent energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (eV)')
    ax.set_xlim(0.0,0.25)
    plt.legend()

    fig.savefig('energy.pdf')

if kinetic:
    array = np.genfromtxt('kinetic.csv', delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1])
    ax.set_title('Three-particle system time-dependent kinetic energy')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Kinetic energy (eV)')
    ax.set_xlim(0.0,0.25)

    fig.savefig('kinetic.pdf')

if kinetic_spectrum:
    array = np.genfromtxt('powerspectrum_kinetic.csv',
    delimiter=',', skip_header=1)

    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(array[:, 0], array[:, 1])
    ax.set_title('Three-particle system kinetic energy spectrum')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Power Spectrum (arb.unit)')
    ax.set_xlim(-100,100)

    fig.savefig('powerspectrum_kinetic.pdf')

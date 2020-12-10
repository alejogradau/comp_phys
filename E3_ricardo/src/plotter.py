#!/usr/bin/env python
###############################################################################
# E3
# Plots a histogram of the generated values for x_i and P(x_i)
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

x_i = np.genfromtxt('./out/generated_xi.csv', delimiter=',', skip_header=1)
P_xi = np.genfromtxt('./out/generated_p_xi.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots(figsize=(10,4))
plt.hist(x_i, bins = 100)
fig.savefig('./out/generated_xi.pdf')

fig, ax = plt.subplots(figsize=(10,4))
plt.hist(P_xi, bins = 100)
fig.savefig('./out/generated_p_xi.pdf')

autocorrelation = np.genfromtxt('./out/autocorrelation.csv', delimiter=',', skip_header=1)
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(autocorrelation[0:400, 0], autocorrelation[0:400, 1])
fig.savefig('./out/autocorrelation.pdf')

autocorrelation = np.genfromtxt('./out/block_averaging.csv', delimiter=',', skip_header=1)
fig, ax = plt.subplots(figsize=(10,4))
ax.plot(autocorrelation[:, 0], autocorrelation[:, 1])
fig.savefig('./out/block_averaging.pdf')

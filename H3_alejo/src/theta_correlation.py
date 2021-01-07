import numpy as np
from numpy import linalg as LA
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

datapath = './out/'
figpath = './plots/'
dataext = '.csv'
figext = '.pdf'

dataname = 'configurations'
figname = 'histogram'

angle = False
cosangle = True

configurations = np.genfromtxt(datapath + dataname + dataext,
                            delimiter=',', skip_header=1)

num_rows, num_cols = configurations.shape

unit_e1 = np.zeros((num_rows, 3))
unit_e2 = np.zeros((num_rows, 3))
x = np.zeros(num_rows)
thetas = np.zeros(num_rows)

e1 = configurations[:, :3]  # Electron1 coordinates
e2 = configurations[:, 3:]  # Electron2 coordinates

norm_e1 = LA.norm(e1, axis=1)
norm_e2 = LA.norm(e2, axis=1)
print(norm_e1.shape)

# Build unit vectors
for i in range(num_rows):
    unit_e1[i, :] = np.divide(e1[i, :], norm_e1[i])
    unit_e2[i, :] = np.divide(e2[i, :], norm_e2[i])

# Build dot product array
for i in range(num_rows):
    x[i] = np.dot(unit_e1[i], unit_e2[i])
print(x)
# Build theta arrays
thetas = np.arccos(x)
print(thetas)
#hist, bin_edges = np.histogram(x, bins=11)
#print(hist, bin_edges)

if cosangle:
    # the histogram of the data
    plt.hist(x, bins=11)
    plt.xlabel(r'$x = \cos(\theta)$')
if angle:
    # the histogram of the data
    plt.hist(thetas, bins=11)
    plt.xlabel(r'$\theta$')

print('Histogram done')


plt.ylabel('Probability')
plt.title('Angular distribution')
plt.savefig(figpath + figname + figext)

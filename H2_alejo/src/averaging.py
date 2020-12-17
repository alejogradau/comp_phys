import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import glob

datapath = './out/task4/np30/'
filelist = glob.glob(datapath+'*.csv')

alpha_mean = 0
alpha2_mean = 0
counter = 0

file = open('alpha_average.txt', "w")

for fname in filelist:
    alpha_trajectory = np.genfromtxt(fname, delimiter=',', skip_header=1)
    alpha = alpha_trajectory[:, 1]

    for i in np.arange(15,29,1):  # Last 15 datapoints
        alpha_mean += alpha[i]
        alpha2_mean += (alpha[i])**2
        counter += 1

alpha_mean /= counter
alpha2_mean /= counter
sigma_alpha = np.sqrt(alpha2_mean - (alpha_mean**2))

header = 'alpha_mean, sigma_alpha\n'

file.write(header)
file.write(str(alpha_mean) + ', ' + str(sigma_alpha))
file.close()

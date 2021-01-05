import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import glob
import decimal

average_E_mean = True
average_alpha = False

if average_E_mean:
    print('We are inside first if')

    def toDecimal(x):  # x is a float.
        return decimal.Decimal(str(x)).quantize(decimal.Decimal('1.00'))
    alphas = np.array([0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23])

    alphas = map(toDecimal, alphas)
    print(alphas)
    header = 'alpha, E_mean_averaged, sigma_E_averaged\n'
    file = open('./out/alpha_vs_E_N1e8.txt', "w")
    file.write(header)
    for i in alphas:  # Calculate E_mean_averaged for such alpha averaging the 5 runs

        counter = 0
        E_mean_averaged = 0
        sigma_E_averaged = 0
        datapath = './out/mean_energy_alpha/mean_energy_alpha' + str(i) + '*'
        print(datapath)
        filelist = glob.glob(datapath)
        #print(filelist)
        for fname in filelist:  # filelist contains 5 files (one per run)
            print(fname)
            data = np.genfromtxt(fname, delimiter=',', skip_header=1)
            E_run = data[0]
            sigma_run = data[1]
            E_mean_averaged += E_run  # Variable to store the average of 5 runs
            sigma_E_averaged += sigma_run
            counter += 1

        E_mean_averaged /= counter
        sigma_E_averaged /= counter
        file.write(str(i) + ', ' + str(E_mean_averaged) + ', ' + str(sigma_E_averaged) + '\n')

    file.close()

if average_alpha:
    alpha_mean = 0
    alpha2_mean = 0
    counter = 0
    file = open('name of file for alpha average.txt', "w")
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

import numpy as np

# Import data
data = np.loadtxt('csv/energy_volume.csv', delimiter=',', skiprows=1)
x = data[:,0]
y = data[:,1]
# print(data, x, y)

# Polynomial fit
coefs = np.polyfit(x, y, 2)

# Create polynomial with coefs
p = np.poly1d(coefs)
xp = np.linspace(64,69,100)

# Estimate the volume of equilibrium
arg = np.argmin(p(xp))
a0 = np.cbrt(xp[arg])
print('The theoretical lattice parameter at T=0K is: ' + str(a0))

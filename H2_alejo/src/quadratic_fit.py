import numpy as np

# Import data
data = np.loadtxt('csv/energy_volume.csv', delimiter=',', skiprows=1)
x = data[:,0]
y = data[:,1]
# print(data, x, y)

# Polynomial fit
coefs, covar = np.polyfit(x, y, 2,  full = False, cov=True)

p0 = coefs[0]
p1 = coefs[1]

err_p0 = covar[0,0]
err_p1 = covar[1,1]

V0 = - (p1)/(2*p0)
err_V0 = V0 * np.sqrt(((err_p0/(2*p0))**2)+((err_p1/(2*p1))**2))

a0 = np.cbrt(V0)
err_a0 = a0 * (err_V0/(3*V0))
print('The theoretical lattice parameter at T=0K is: ' + str(a0)
+ ' +- ' + str(err_a0))

# Create polynomial with coefs
#p = np.poly1d(coefs)
#xp = np.linspace(64,69,100)

# Estimate the volume of equilibrium
#arg = np.argmin(p(xp))
#a0 = np.cbrt(xp[arg])
#print('The theoretical lattice parameter at T=0K is: ' + str(a0))

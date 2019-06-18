#!/usr/bin/python
import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from scipy.constants import e,k, pi, epsilon_0
import math

# mass number
A = 27
# decrease in KE per collision
a = ((1-A)/(1+A))**2
# atomic number
Z = 13
# coulomb constant
Q = 1/(4*pi*epsilon_0)
# Temperature of the plasma
T = 1000
# number density of aluminium (FCC)
N = 6 * 10 **28
# density of aluminium
rho = 2710
# di inner thickness
di = 1 * 10**(-6)
U = 3*k*T * 10 **5 / 2
C = math.sqrt(6*k*T)/(20*N*pi*Z**2 * Q**2 * e **4)
B = 2*pi*e*U/(rho * di)



energy = np.arange(100, 100000,10)
#n = np.zeros(shape = energy.shape, dtype = int)
#l = np.zeros(shape = energy.shape, dtype = float)
# calculate thermal energy
TE = k * 500/ e

def collision(E):
    if E <= TE:
        return 0
    else:
        return math.ceil(math.log(TE / E, a))

#def path(n, E):
#    return (1-a**(3*(n+1)/2))/(1-a**(3/2))*E**(3/2)

def integrand2(E):
    if E <= TE:
        n = 0
    else:
        n = math.ceil(math.log(TE / E, a))
    l = C*(1-a**(3*(n+1)/2))/(1-a**(3/2))*(e*E - e*U)**(3/2)
    return B*((math.sqrt((e*U/e*E) + 4*di**2 * (1-e*U/e*E)/l**2) - math.sqrt(e*U/e*E)))*10 **6 * math.exp(-E/120000)

def integrand1(E):
    return B*(1-math.sqrt(e*U/e*E))* 10 **6

I = np.array([integrand2(E) for E in energy])
integral2 = quad(integrand2, 1000000, np.inf)
integral = quad(integrand1, e*U, 1000000)
integral2 = [x/e for x in integral2]
integral = [x/e for x in integral]
print(tuple(map(sum, zip(integral, integral2))))

#n = np.array([collision(E) for E in energy])
'''
for i in range(n.shape[0]):
    l[i] = (1-a**(3*(n[i]+1)/2))/(1-a**(3/2))*energy[i]**(3/2)

integrand = np.zeros(shape = energy.shape, dtype=float)

for i in range(n.shape[0]):
    integrand[i] = - math.sqrt(1/energy[i]) * math.exp(-energy[i])
    math.sqrt((1/energy[i]) + 4*e*(1-e/energy[i])/(l[i]**2))

'''
plt.plot(energy, I)
plt.show()

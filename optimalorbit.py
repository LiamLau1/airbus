#!/usr/bin/python
from ctypes import *
from math import acos
from numpy import sqrt, sin, cos, pi
import scipy.integrate as integrate
from mpl_toolkits.mplot3d import Axes3D


# Loading the standard c library
cdll.LoadLibrary("libc.so.6")
libc = CDLL("libc.so.6")
# We first load the precompiled AP8 Model C shared library into python so we can call the function. Here we use the data for AP8MAX
simlib = CDLL("/home/liamlhlau/Documents/airbus/simulation/libmodel.so")
# We need to change the return type for ctypes
simlib.simulate.restype = c_double

# Orbit code
def makecubelimits(axis, centers=None, hw=None):
    lims = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    if centers == None:
        centers = [0.5*sum(pair) for pair in lims]

    if hw == None:
        widths  = [pair[1] - pair[0] for pair in lims]
        hw      = 0.5*max(widths)
        ax.set_xlim(centers[0]-hw, centers[0]+hw)
        ax.set_ylim(centers[1]-hw, centers[1]+hw)
        ax.set_zlim(centers[2]-hw, centers[2]+hw)
        #print("hw was None so set to:", hw)
    else:
        try:
            hwx, hwy, hwz = hw
            print("ok hw requested: ", hwx, hwy, hwz)

            ax.set_xlim(centers[0]-hwx, centers[0]+hwx)
            ax.set_ylim(centers[1]-hwy, centers[1]+hwy)
            ax.set_zlim(centers[2]-hwz, centers[2]+hwz)
        except:
            print("nope hw requested: ", hw)
            ax.set_xlim(centers[0]-hw, centers[0]+hw)
            ax.set_ylim(centers[1]-hw, centers[1]+hw)
            ax.set_zlim(centers[2]-hw, centers[2]+hw)

    return centers, hw

TLE = """1 43205U 18017A   18038.05572532 +.00020608 -51169-6 +11058-3 0  9993
2 43205 029.0165 287.1006 3403068 180.4827 179.1544 08.75117793000017"""
L1, L2 = TLE.splitlines()

from skyfield.api import Loader, EarthSatellite
from skyfield.timelib import Time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

halfpi, pi, twopi = [f*np.pi for f in [0.5, 1, 2]]
degs, rads = 180/pi, pi/180

load = Loader('~/Documents/fishing/SkyData')
data = load('de421.bsp')
ts   = load.timescale()

planets = load('de421.bsp')
earth   = planets['earth']

Roadster = EarthSatellite(L1, L2)

#print(Roadster.epoch.tt)
hours = np.arange(0, 3, 0.01)

time = ts.utc(2018, 2, 7, hours)

Rpos    = Roadster.at(time).position.km
Rposecl = Roadster.at(time).ecliptic_position().km

#print(Rpos.shape)

re = 6378.

theta = np.linspace(0, twopi, 201)

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
#U = [3*k*T * f/ 2 for f in flux]
C = 2* math.sqrt(6*k*T)/(5*N*pi*Z**2 * Q**2 * e **4)
#B = [2*pi*e*pot/(rho * di) for pot in U] #Find all the constant factors based off flux
#m = 10 **4 #plasma density parameter (nplasma density divided by thickness of shell)



#energy = np.arange(100, 100000,10)
#n = np.zeros(shape = energy.shape, dtype = int)
#l = np.zeros(shape = energy.shape, dtype = float)
# calculate thermal energy
TE = k * 500/ e

def collision(E):
    return 0 if E <= TE else math.ceil(math.log(TE / E, a))

#def path(n, E):
#    return (1-a**(3*(n+1)/2))/(1-a**(3/2))*E**(3/2)

def integrand2(E, f, p, b):
    if E <= TE:
        n = 0
    else:
        n = math.ceil(math.log(TE / E, a))
    l = C*(1-a**(3*(n+1)/2))/(1-a**(3/2))*(e*E - e*p)**(3/2)
    return b*((math.sqrt((e*p/e*E) + 4*di**2 * (1-e*p/e*E)/l**2) - math.sqrt(e*p/e*E)))*(f/120000) * math.exp(-E/120000) # don't need e * as we're integrating dE

def integrand1(E, f, p, b):
    return b*(1-math.sqrt(e*p/e*E))*(f/120000) * math.exp(-E/120000)  # don't need e* as we're integrating dE, 120000

def g(r):
    # we now numerically integrate with differing flux values due to changing r We now have args to show the changing parameter
    f = simlib.simulate(c_double(r),c_double(0)) # due to the axisymmetry we have lattitude 0
    U = 3*k*T * f /(2)  #U as a function of f 
    B = 2*pi*e*U/(rho * di) #constant as a function of the parameter U
    integral2 = quad(integrand2, 100000, np.inf, args=(f,U,B)) # I'm not sure what the root E* is, initially say 100000
    integral = quad(integrand1, e*U, 100000, args=(f,U,B)) # I'm not sure what the root E* is
    integral2 = [x/e for x in integral2] # converting into watts/kg
    integral = [x/e for x in integral] # converting into watts/kg
    w = integral[0] + integral2[0]
    return w

#print(tuple(map(sum, zip(integral, integral2))))
cth, sth, zth = [f(theta) for f in [np.cos, np.sin, np.zeros_like]]
lon0 = re*np.vstack((cth, zth, sth))
lons = []
for phi in rads*np.arange(0, 180, 15):
    cph, sph = [f(phi) for f in [np.cos, np.sin]]
    lon = np.vstack((lon0[0]*cph - lon0[1]*sph,
                     lon0[1]*cph + lon0[0]*sph,
                     lon0[2]) )
    lons.append(lon)

lat0 = re*np.vstack((cth, sth, zth))
lats = []
for phi in rads*np.arange(-75, 90, 15):
    cph, sph = [f(phi) for f in [np.cos, np.sin]]
    lat = re*np.vstack((cth*cph, sth*cph, zth+sph))
    lats.append(lat)

#fig = plt.figure(figsize=[10, 8])  # [12, 10]

#ax  = fig.add_subplot(1, 1, 1, projection='3d')


r_Roadster = np.sqrt((Rpos**2).sum(axis=0))
alt_roadster = r_Roadster - re

### Main algorithm ###
#from the specific power orbit we know the peak is between this range
up_lim = 30000
low_lim = 10000
iterationN = 4 # number of iterations- converges to 3dp in 3 iterations

# take second element for sort
def takeSecond(elem):
    return elem[1]

def iterate(n,a,b):
    if n == 0:
        return a
    else:
        rs = np.linspace(a,b,1000)
        ws = list(map(g, rs)) # compute all the ws for all the rs, need list to convert map and zip objects into lists
        pairs = list(zip(rs,ws)) #create tuples of rs and ws
        rstemp = [i[0] for i in sorted(pairs, key=takeSecond, reverse = True)[0:2]] #sorts (rs,ws) by high to low ws then taking the first two tuples and then placing their rs into a list
        r0 = rstemp[0]
        r1 = rstemp[1]
        print(n)
        return iterate(n-1,r0,r1)

def optimal(a,b):
    rs = np.linspace(a,b,1000)
    ws = list(map(g, rs)) # compute all the ws for all the rs, need list to convert map and zip objects into lists
    pairs = list(zip(rs,ws)) #create tuples of rs and ws
    rstemp = [i[0] for i in sorted(pairs, key=takeSecond, reverse = True)[0:2]] #sorts (rs,ws) by high to low ws then taking the first two tuples and then placing their rs into a list
    r0 = rstemp[0]
    r1 = rstemp[1]
    #print(n)
    return r0,r1

r0new = low_lim
r1new = up_lim
while iterationN != 0:
    r0new, r1new = optimal(r0new,r1new)
    iterationN -= 1
###########################################
if iterationN== 0:
    print((r0new+r1new)/2)


x, y, z = Rpos
optimalr = (r0new+r1new)/2
x = np.linspace(-optimalr,optimalr,1000)
xtemp = np.linspace(-optimalr,optimalr,1000)
y = [sqrt((optimalr**2) - xpos**2 )for xpos in x]
ytemp = [-sqrt((optimalr**2) - xpos**2 )for xpos in x]
x = np.concatenate((x,xtemp),axis=0)
y = np.concatenate((y,ytemp),axis=0)
z = np.zeros(x.shape[0], dtype=float)
radius = []
phi_lat = []

for i in range(0,x.shape[0]):#changed from Rpos.shape[1], x.shape[0]
    radii = sqrt(x[i]**2 + y[i]**2 + z[i]**2)
    radius.append(radii)
for i in range(0,x.shape[0]):
    lat_calc = pi/2 - acos(z[i]/radius[i])
    phi_lat.append(lat_calc)

flux = []
Lvalue = [] # McIlwain Parameter

for i in range(0,x.shape[0]):
    altrad = phi_lat[i]*pi/180.0 #magnetic latitude, not latitude ( +/- 11 deg )
    # get the field by simple dipole, value in nanoteslas
    BetaValueF = 30610.0*sqrt(cos(1.57-altrad)*cos(1.57-altrad)*3.0+1.0)/(radius[i]**3)
    #compute L value
    bottomF = radius[i]**6
    bottomF = 4 - (BetaValueF*BetaValueF*bottomF/(3.06e4*3.06e4))
    Lvalue.append(3*sqrt(radius[i]**2)/(bottomF * 10**4)) #added 10^4 to normalise to within the accepted values of L
for r in radius:
    flux.append(simlib.simulate(c_double(r),c_double(phi_lat[i])))

w = np.array([])
wvalue = g(optimalr)
print(wvalue)
for k in radius:
    w = np.append(w,wvalue)

#print(iterate(iterationN, low_lim,up_lim))






fig = plt.figure(figsize=[10, 8])  # [12, 10]
ax  = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot(x, y, z)
img = ax.scatter(x, y, z, c=w)
for x, y, z in lons:
    ax.plot(x, y, z, '-k')
for x, y, z in lats:
    ax.plot(x, y, z, '-k')
fig.colorbar(img)

plt.show()



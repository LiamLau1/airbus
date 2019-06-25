#!/usr/bin/python
from ctypes import *
from math import acos, atan2
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
if __name__=="__main__":
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

    if True:
        #fig = plt.figure(figsize=[10, 8])  # [12, 10]

        #ax  = fig.add_subplot(1, 1, 1, projection='3d')

        #x, y, z = Rpos
        sa = 15900
        sb = 14240
        sc = sqrt((sa**2) - (sb**2))
        x = np.concatenate((np.linspace(-sa,lon.min() + sc, 1000),np.linspace(lon.max() - sc,sa, 1000)),axis=0)
        #x = x[(abs(x) <  8500)]
        #x = np.delete(x,[0,1])
        #xtemp = np.delete(xtemp,0)
        z = [+sqrt((sb**2) - (sb**2/sa **2) *(xpos ** 2)) for xpos in x]
        ztemp = [-k for k in z]
        z = np.asarray(z)
        ztemp = np.asarray(ztemp)
        #x = np.concatenate((x,xtemp), axis = 0)
        x = [pk - sc for pk in x]
        xtemp = x
        x = np.asarray(x)
        xtemp = np.asarray(xtemp)
        #xsk = np.concatenate((x,xtemp),axis = 0)
        #zsk = np.concatenate((z,ztemp),axis = 0)
        #z = np.concatenate((z,ztemp), axis = 0)
        y = np.zeros(x.shape[0], dtype=float)
        theta1 = [atan2(pz,px) for px,pz in zip(x,z)]
        theta1 = np.asarray(theta1)
        theta2 = [atan2(pz2,px2) for px2,pz2 in zip(xtemp,ztemp)]
        theta2 = np.asarray(theta2)


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
            #Lvalue.append(3*sqrt(radius[i]**2)/(bottomF * re)) #added 10^4 to normalise to within the accepted values of L
            Lvalue.append(radius[i]/(re *(cos(altrad)**2)))
        for r, phi in zip(radius, phi_lat):
            flux.append(simlib.simulate(c_double(r),c_double(phi))) #for earth
        flux = np.asarray(flux)
        Lvalue = np.asarray(Lvalue)

        radiustemp = []
        phi_lattemp = []
        for i in range(0,xtemp.shape[0]):#changed from Rpos.shape[1], x.shape[0]
            radiitemp = sqrt(xtemp[i]**2 + y[i]**2 + ztemp[i]**2)
            radiustemp.append(radiitemp)
        for i in range(0,xtemp.shape[0]):
            lat_calctemp = pi/2 - acos(ztemp[i]/radius[i]) 
            phi_lattemp.append(lat_calctemp)

        fluxtemp = []
        Lvaluetemp = [] # McIlwain Parameter

        for i in range(0,xtemp.shape[0]):
            altradtemp = phi_lattemp[i]*pi/180.0 #magnetic latitude, not latitude ( +/- 11 deg ) 
            # get the field by simple dipole, value in nanoteslas
            BetaValueFtemp = 30610.0*sqrt(cos(1.57-altradtemp)*cos(1.57-altradtemp)*3.0+1.0)/(radiustemp[i]**3)
            #compute L value 
            bottomFtemp = radiustemp[i]**6
            bottomFtemp = 4 - (BetaValueFtemp*BetaValueFtemp*bottomFtemp/(3.06e4*3.06e4))
            #Lvalue.append(3*sqrt(radius[i]**2)/(bottomF * re)) #added 10^4 to normalise to within the accepted values of L
            Lvaluetemp.append(radiustemp[i]/(re *(cos(altradtemp)**2)))
        for r, phi in zip(radiustemp, phi_lattemp):
            fluxtemp.append(simlib.simulate(c_double(r),c_double(phi))) #for earth
        fluxtemp = np.asarray(fluxtemp)
        Lvaluetemp = np.asarray(Lvaluetemp)




        '''
        ax.plot(x, y, z)
        for x, y, z in lons:
            ax.plot(x, y, z, '-k')
        for x, y, z in lats:
            ax.plot(x, y, z, '-k')
        '''

        #centers, hw = makecubelimits(ax)

        #print("centers are: ", centers)
        #print("hw is:       ", hw)

        #plt.show()

    r_Roadster = np.sqrt((Rpos**2).sum(axis=0))
    alt_roadster = r_Roadster - re


    ''' 
     if True:
        plt.figure()
        plt.plot(hours, r_Roadster)
        plt.plot(hours, alt_roadster)
        plt.xlabel('hours', fontsize=14)
        plt.ylabel('Geocenter radius or altitude (km)', fontsize=14)
        plt.show()
    '''

import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
from scipy.constants import e,k, pi, epsilon_0, G
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
T = 5000 #PARAMETER
# number density of aluminium (FCC)
N = 6 * 10 **28
# density of aluminium
rho = 2710
# di inner thickness
di = 1 * 10**(-6)
#U = [3*k*T * f/ 2 for f in flux]
C = 2 * math.sqrt(6*k*T)/(5*N*pi*Z**2 * Q**2 * e **4)
#B = [2*pi*e*pot/(rho * di) for pot in U] #Find all the constant factors based off flux
#m = 10 **2 #correction



#energy = np.arange(100, 100000,10)
#n = np.zeros(shape = energy.shape, dtype = int)
#l = np.zeros(shape = energy.shape, dtype = float)
# calculate thermal energy
TE = k * 500/ e

def collision(E):
    return 0 if E <= TE else math.ceil(math.log(TE / E, a))
#def path(n, E):
#    return (1-a**(3*(n+1)/2))/(1-a**(3/2))*E**(3/2)

def integrand2(E, f, p, b,E0):
    if E <= TE:
        n = 0
    else:
        n = math.ceil(math.log(TE / E, a))
    l = C*(1-a**(3*(n+1)/2))/(1-a**(3/2))*(e*E - e*p)**(3/2)
    return b*((math.sqrt((e*p/e*E) + 4*di**2 * (1-e*p/e*E)/l**2) - math.sqrt(e*p/e*E)))*(f/E0) * math.exp(-E/E0) # don't need e * as we're integrating dE

def integrand1(E, f, p, b,E0):
    return b*(1-math.sqrt(e*p/e*E))*(f/E0) * math.exp(-E/E0)  # don't need e* as we're integrating dE, E0

w = np.array([])
# we now iterate at each location the integral with differing flux values which changes U, B and thus the integrands. We now have args to show the changing parameter
for f,L in zip(flux,Lvalue):
    U = 3*k*T * f /(2)  #U as a function of f
    B = 2*pi*e*U/(rho * di) #constant as a function of the parameter U
    M = 5**3 * 120000 #constant for E0
    E0temp = M * (L**-3)
    integral2 = quad(integrand2, 100000, np.inf, args=(f,U,B,E0temp)) # I'm not sure what the root E* is, initially say 100000
    integral = quad(integrand1, e*U, 100000, args=(f,U,B,E0temp)) # I'm not sure what the root E* is
    integral2 = [x/e for x in integral2] # converting into watts/kg
    integral = [x/e for x in integral] # converting into watts/kg
    w = np.append(w, integral[0] + integral2[0])
wtemp = np.array([])
# we now iterate at each location the integral with differing flux values which changes U, B and thus the integrands. We now have args to show the changing parameter
for f,L in zip(fluxtemp,Lvaluetemp):
    U = 3*k*T * f /(2)  #U as a function of f
    B = 2*pi*e*U/(rho * di) #constant as a function of the parameter U
    M = 5**3 * 120000 #constant for E0
    E0temp = M * (L**-3)
    integral2temp = quad(integrand2, 100000, np.inf, args=(f,U,B,E0temp)) # I'm not sure what the root E* is, initially say 100000
    integraltemp = quad(integrand1, e*U, 100000, args=(f,U,B,E0temp)) # I'm not sure what the root E* is
    integral2temp = [x/e for x in integral2temp] # converting into watts/kg
    integraltemp = [x/e for x in integraltemp] # converting into watts/kg
    wtemp = np.append(wtemp, integraltemp[0] + integral2temp[0])

ME = 5.972 * 10**24
TimeP = 2 * pi * sqrt(((sa*10**3) **3)/(G * ME))
#omega = 2* pi /TimeP
ecc = math.sqrt(1-((sb**2)/(sa**2)))
Eori = [acos((1-(r/sa))/ecc) for r in radius]
Etemp = [acos((1-(r2/sa))/ecc) for r2 in radiustemp]
tori = [TimeP/(2*pi)*(E-ecc*cos(E)) for E in Eori]
ttemp = [TimeP/(2*pi)*(E-ecc*cos(E)) for E in Etemp]
#print(tuple(map(sum, zip(integral, integral2))))

#### Code to make equal aspect ratio
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz)) * 6 * 10**4
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


print(flux.shape[0])
print(Lvalue.shape[0])
print(w.shape[0])


#fig = plt.figure(figsize=[10, 8])  # [12, 10]
#ax  = fig.add_subplot(1, 1, 1, projection='3d')
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.set_aspect(aspect='equal')
axisEqual3D(ax)
ax.set_xlabel('x/ km', fontsize = 15)
ax.set_ylabel('y/ km', fontsize = 15)
ax.set_zlabel('z/ km', fontsize = 15)
ax.plot(x, y, z)
ax.plot(xtemp,y,ztemp)
img = ax.scatter(x, y, z, c=w)
img2 = ax.scatter(xtemp,y,ztemp, c = wtemp)


for x, y, z in lons:
    ax.plot(x, y, z, '-k')
for x, y, z in lats:
    ax.plot(x, y, z, '-k')
fig.colorbar(img)
fig.colorbar(img2)


#plt.figure()
#plt.plot(Lvalue,flux, '.')
#plt.plot(radius,flux, '.')
plt.figure()
plt.plot(theta1,w,'b')
plt.plot(theta2,wtemp,'b')
plt.figure()
plt.plot(tori,w,'r')
plt.figure()
plt.plot(ttemp,wtemp,'r')
plt.show()


print((w.mean(), w.std()/w.shape[0]))
print(w.max())
Lvalue = np.asarray(Lvalue)
print((Lvalue.mean(), Lvalue.std()/Lvalue.shape[0]))


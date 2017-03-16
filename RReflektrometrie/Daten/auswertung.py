# -*- coding: utf-8 -*-

import numpy as np
from numpy.lib import scimath
import scipy
from scipy.optimize import curve_fit
from scipy import constants
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import log10,log,exp,sqrt
import sympy
import cmath

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

data_signal = np.genfromtxt('messung1_1306.uxd', unpack='True', skip_header=66)

data_noise = np.genfromtxt('messung2_1349.uxd', unpack='True', skip_header=66)

data = np.array([data_signal[0],data_signal[1]-data_noise[1]])

datatxt = open('messung_bereinigt.txt', 'w')

for i in range (0,len(data[0,:])):

    datatxt.write(str(data[0,i]) + '\t' + str(data[1,i]) + '\n')

datatxt.close()

plt.figure()
plt.plot(data_signal[0], data_signal[1], 'r-', label="Messdaten")
plt.plot(data_noise[0], data_noise[1], 'g-', label="diffuser Scan")
plt.plot(data[0], data[1], 'b-', label="bereinigte Messdaten")
plt.yscale("log")
plt.xlabel(r"$\alpha_{i} / ^{\circ}$")
plt.ylabel("Intensität / a.u.")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("../Protokoll/images/rawdata.pdf")


l=1.54*10**(-10) #wavelength in Angstrom
k=2*np.pi/l #wavenumber
ai = np.linspace(0,2.5,10000) #incident angle
qz = 4*np.pi/(l*10**(10))*np.sin(np.deg2rad(ai)) #change of wave vector
n1 = 1 #n(air)
ag = 0.3247 #geometry-angle

#Peak-Detection

data_schwing_mod = data[:,50:312]

peaks = np.ones((2,26))

for i in range (1,26):

     m = (i-1)*10
     n = m+10
     index = np.argmin(data_schwing_mod[1,m:n])
     index=index+m
     peaks[0,i] = data_schwing_mod[0,index]
     peaks[1,i] = data_schwing_mod[1,index]

peaks = np.delete(peaks, 0, 1)
peaks = np.delete(peaks, 10, 1)
dist = np.ones(23)

for i in range (1,23):

    dist[i]=peaks[0,i+1]-peaks[0,i]

dist = np.delete(dist, [0,15], 0)

mean_dist = ufloat(np.mean(np.deg2rad(dist)),np.std(np.deg2rad(dist)))
z_unp = l/(2*mean_dist)
z2 = unp.nominal_values(z_unp)

print("{:.2u}".format(z_unp))
print(dist)

def geometryfactor(I,a):

    if (a<=ag):
        G=np.sin(np.deg2rad(a))/np.sin(np.deg2rad(ag))
    else:
        G=1

    return G*I

for i in range (0,len(data[0])):
    data[1,i]=geometryfactor(data[1,i],data[0,i])


def reflectrometry(a,n2,n3,s1,s2):

    a0=5*10**7

    kz1 =k*scimath.sqrt(n1**2-(np.cos(a))**2)
    kz2 =k*scimath.sqrt(n2**2-(np.cos(a))**2)
    kz3 =k*scimath.sqrt(n3**2-(np.cos(a))**2)

    r12 = (kz1-kz2)/(kz1+kz2)*np.exp(-2*kz1*kz2*s1**2)
    r23 = (kz2-kz3)/(kz2+kz3)*np.exp(-2*kz2*kz3*s2**2)
    x2 = np.exp(-2*1j*kz2*0.96*z2)*r23
    x1 = (r12+x2)/(1+r12*x2)
    rr = a0*(np.absolute(x1))**2
    return rr

def electron_density(de):

    return (2*np.pi)/(l**2*constants.value(u'classical electron radius'))*de


x0 = [1-25*10**(-7), 1-70*10**(-7), 75*10**(-11), 35*10**(-11)]

params, covariance = curve_fit(reflectrometry, np.deg2rad(data[0,40:300]), data[1,40:300], p0=x0)
errors = np.sqrt(np.diag(covariance))

for i in range (0,len(params)):
    print(params[i], "\t", "+/-", errors[i], "\n")

delta = unp.uarray([params[0],params[1]],[errors[0], errors[1]])
#calculation electron-density

print('Die Elektronendichten lauten: ', '\n')
print(electron_density(1-delta[0]), '\n')
print(electron_density(1-delta[1]), '\n')

plt.figure()
plt.plot(data[0], data[1], 'b-', label="Messwerte")
plt.plot(peaks[0], peaks[1], 'rx', label="verwendete Minima")
plt.plot(ai, reflectrometry(np.deg2rad(ai), *params), 'y-', label="Fit")
#plt.plot(ai, reflectrometry(np.deg2rad(ai), *x0), 'g-', label="Test")
plt.yscale("log")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("../Protokoll/images/reflectrometry.pdf")

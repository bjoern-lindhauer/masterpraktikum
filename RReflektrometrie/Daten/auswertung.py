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


l=1.54*10**(-10) #wavelength in Angstrom
k=2*np.pi/l #wavenumber
ai = np.linspace(0,2.5,10000) #incident angle
qz = 4*np.pi/(l*10**(10))*np.sin(np.deg2rad(ai)) #change of wave vector
n1 = 1 #n(air)

#Peak-Detection

data_schwing_mod = data[:,50:312]

peaks = np.ones((2,26))

for i in range (1,26):

     m = (i-1)*10
     n = m+10
     index = np.argmax(data_schwing_mod[1,m:n])
     index=index+m
     peaks[0,i] = data_schwing_mod[0,index]
     peaks[1,i] = data_schwing_mod[1,index]

peaks = np.delete(peaks, 0, 1)
peaks = np.delete(peaks, 16, 1)
dist = np.ones(23)

for i in range (1,23):

    dist[i]=peaks[0,i+1]-peaks[0,i]

dist = np.delete(dist, [0,15], 0)

z2 = l/(2*np.mean(np.deg2rad(dist)))

print(z2)

def reflectrometry(a,n2,n3,s1,s2,a0):

    kz1 =k*scimath.sqrt(n1**2-(np.cos(a))**2)
    kz2 =k*scimath.sqrt(n2**2-(np.cos(a))**2)
    kz3 =k*scimath.sqrt(n3**2-(np.cos(a))**2)

    r12 = (kz1-kz2)/(kz1+kz2)*np.exp(-2*kz1*kz2*s1**2)
    r23 = (kz2-kz3)/(kz2+kz3)*np.exp(-2*kz2*kz3*s2**2)
    x2 = np.exp(-2*1j*kz2*z2)*r23
    x1 = (r12+x2)/(1+r12*x2)
    rr = a0*(np.absolute(x1))**2
    return rr


x0 = [1-20*10**(-7), 1-75*10**(-7), 55*10**(-11), 35*10**(-11), 3*10**7]

params, covariance = curve_fit(reflectrometry, np.deg2rad(data[0,0:300]), data[1,0:300], p0=x0)
errors = np.sqrt(np.diag(covariance))

for i in range (0,len(params)):
    print(params[i], "\t", "+/-", errors[i], "\n")


plt.figure()
plt.plot(data[0], data[1], 'b-', label="Messwerte")
plt.plot(peaks[0], peaks[1], 'rx', label="verwendete Peaks")
#plt.plot(ai, reflectrometry(np.deg2rad(ai), *params), 'y-', label="Fit")
plt.plot(ai, reflectrometry(np.deg2rad(ai), *x0), 'g-', label="Test")
plt.yscale("log")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("reflectrometry.pdf")

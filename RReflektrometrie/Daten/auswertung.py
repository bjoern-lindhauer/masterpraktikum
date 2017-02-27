# -*- coding: utf-8 -*-

import numpy as np
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

ai = np.linspace(0,5,10000)*np.pi/180 #incident angle
qz = 4*np.pi/(l*10**(10))*np.sin(ai) #change of wave vector


def reflectrometry(a,n1,n2,n3,s1,s2,a0):

    #z-components
    kz1=k*np.sqrt(n1**2-np.square(np.cos(a)))
    kz2=k*np.sqrt(n2**2-np.square(np.cos(a)))
    kz3=k*np.sqrt(n3**2-np.square(np.cos(a)))

    #modified Fresnel-coefficients

    r12=(kz1-kz2)/(kz1+kz2)*np.exp(-2*kz1*kz2*s1**2)
    r23=(kz2-kz3)/(kz2+kz3)*np.exp(-2*kz2*kz3*s2**2)
    x2=np.exp(-2*1j*kz2*z2)*r23
    x1=a0*(r12+x2)/(1+r12*x2)

    return np.absolute(x1)

#Abfallende Schwingung

data_schwing_mod = data[:,50:312]

#Peak-Detection

peaks = np.ones((2,26))

for i in range (1,26):

     m = (i-1)*10
     k = m+10
     index = np.argmax(data_schwing_mod[1,m:k])
     index=index+m
     peaks[0,i] = data_schwing_mod[0,index]
     peaks[1,i] = data_schwing_mod[1,index]

peaks = np.delete(peaks, 0, 1)

dist = np.ones(24)

for i in range (1,24):

    dist[i]=peaks[0,i+1]-peaks[0,i]

dist = np.delete(dist, [0,15], 0)

z2 = l/(2*np.mean(dist))

x0 = [1, 1-10**(-6), 1-10**(-6), 0.5*10**(-9), 0.5*10**(-9), 7*10**8]

params, covariance = curve_fit(reflectrometry, data[0,:], data[1,:], p0=x0)



plt.figure()
plt.plot(data[0], data[1], 'b-', label="Messwerte")
plt.plot(ai, reflectrometry(ai, *params), 'g-', label="Fit")
plt.yscale("log")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("reflectrometry.pdf")

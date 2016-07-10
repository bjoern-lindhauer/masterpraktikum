# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 15:07:54 2015

@author: robert
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy import constants
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import log10,log,exp,sqrt
import sympy

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

#Funktionen definieren

def f(x,m, b):
    return m*x+b

def magnhelm(I, N, R):
    return constants.mu_0*(8*I*N)/(np.sqrt(125)*R)

#Daten auslesen

daten1 = np.genfromtxt('landebestimmung.txt', unpack='True')
daten2 = np.genfromtxt('periodenauswertung.txt', unpack='True')

#Daten transformieren

B1 = magnhelm(daten1[1]*0.1, 11, 0.1639) + magnhelm(daten1[3]*0.3, 154, 0.1579)
B2 = magnhelm(daten1[2]*0.1, 11, 0.1639) + magnhelm(daten1[4]*0.3, 154, 0.1579)

#Erste Magnetfeld-Messung fitten

x_plot=np.linspace(100,1000, num=10000)
params, covariance = curve_fit(f,daten1[0], B1)

#Daten und Fit plotten

plt.figure()

plt.plot(daten1[0], B1, 'bx', label="Messung erstes Isotop")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='linearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(100,1010)
#plt.ylim(0,6)
plt.xlabel(r'Frequenz $\nu$ [kHz]')
plt.ylabel(r'Magnetfeld B [T]')
plt.savefig('plot1a.pdf')

errors=np.sqrt(np.diag(covariance))

print('m =', params[0], '+/-', errors[0])
print('b =', params[1], '+/-', errors[1])

# Zweite Magnetfeld-Messung fitten

params, covariance = curve_fit(f,daten1[0], B2)

#Daten und Fit plotten

plt.figure()

plt.plot(daten1[0], B2, 'bx', label="Messung zweites Isotop")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='linearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(100,1010)
#plt.ylim(0,6)
plt.xlabel(r'Frequenz $\nu$ [kHz]')
plt.ylabel(r'Magnetfeld B [T]')
plt.savefig('plot1b.pdf')

errors=np.sqrt(np.diag(covariance))

print('m =', params[0], '+/-', errors[0])
print('b =', params[1], '+/-', errors[1])

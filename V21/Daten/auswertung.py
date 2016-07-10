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

def gi(m):
    return (4*np.pi*m0)/(e0*m)*10**3                    #10**3 wg. kHz#

def I(gj):
    return (g/(4*gj))-1+sqrt((1-(g/(4*gj)))**2-(3/4)*(1-(g/gj)))


#Daten auslesen

daten1 = np.genfromtxt('landebestimmung.txt', unpack='True')
daten2 = np.genfromtxt('periodenauswertung.txt', unpack='True')

#Daten transformieren

B1 = magnhelm(daten1[1]*0.1, 11, 0.1639) + magnhelm(daten1[3]*0.3, 154, 0.1579)
B2 = magnhelm(daten1[2]*0.1, 11, 0.1639) + magnhelm(daten1[4]*0.3, 154, 0.1579)

#Magnetfeld-Messung fitten

x_plot=np.linspace(100,1000, num=10000)
params, covariance = curve_fit(f,daten1[0], B1)
params2, covariance2 = curve_fit(f,daten1[0], B2)

#Daten und Fit plotten

plt.figure()

plt.plot(daten1[0], B1, 'bx', label="Messung erstes Isotop")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='linearer Fit erstes Isotop')
plt.plot(daten1[0], B2, 'kx', label="Messung zweites Isotop")
plt.plot(x_plot, f(x_plot, *params2), 'g-', label='linearer Fit zweites Isotop')
plt.legend(loc="best", numpoints=1)
plt.xlim(100,1010)
#plt.ylim(0,6)
plt.xlabel(r'Frequenz $\nu$ [kHz]')
plt.ylabel(r'Magnetfeld B [T]')
plt.savefig('plot1.pdf')

errors=np.sqrt(np.diag(covariance))

# print('m =', params[0], '+/-', errors[0])
# print('b =', params[1], '+/-', errors[1])

errors2=np.sqrt(np.diag(covariance2))

# print('m =', params2[0], '+/-', errors2[0])
# print('b =', params2[1], '+/-', errors2[1])

#gi berechnen

m1 = ufloat(params[0], errors[0])
m2 = ufloat(params2[0], errors2[0])

m0=9.109*10**(-31)
e0 = constants.elementary_charge

g1 = gi(m1)
g2 = gi(m2)

#Kernspin berechnen

g = 2.002

I1 = I(g1)
I2 = I(g2)

print(I1, I2)
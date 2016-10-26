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

#Daten auslesen

wave = np.genfromtxt('wavelength.txt', unpack='True')

polarisation = np.genfromtxt('polarisation.txt', unpack='True')

tem00 = np.genfromtxt('tem00.txt', unpack='True')

#Fehlerformelausgabe

def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()

    if err_vars == None:
        err_vars = f.free_symbols

    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'

    return latex(sympy.sqrt(s), symbol_names=latex_names)

#Funktionen definieren

def I(r, I0, r0, w):
    return I0*np.exp(-2(r-r0)/w**2)

def Lambda(d,L,n,g):
    return (np.sin(np.arctan(d/L)))/(n*g)

def Ip(phi,I0,delta):
    return I0*np.cos(phi+delta)**2

#Funktion fitten

guess = [1, 1]
params_00, covariance_00 = curve_fit(I, tem00)

#Daten und Fit plotten

#x_plot=np.linspace(0,10, num=1000)
plt.figure()
# err1 = unp.std_devs(data1)
# err2 = unp.std_devs(data2)
#
#
# plt.errorbar(unp.nominal_values(data1) + err1, -unp.nominal_values(data2) + err2, fmt='bx', label="Messung")
plt.plot(x_plot, f(x_plot, 1.4), 'r-', label='Planarer und r=1m Spiegel')
plt.plot(x_plot, g(x_plot, 1, 1.4), 'gx', label='r=1m Spiegel und r=1.4m Spiegel')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,2)
plt.ylim(0,1)
plt.xlabel(r'L/m')
plt.ylabel(r'g$_1*$g$_2$')
plt.show()

# errors=np.sqrt(np.diag(covariance))
#
# print('m =', params[0], '+/-', errors[0])
# print('b =', params[1], '+/-', errors[1])
#

#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

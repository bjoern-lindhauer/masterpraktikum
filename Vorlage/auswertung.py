# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 15:07:54 2015

@author: robert
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import log10,log,exp,sqrt
import sympy

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

#Daten auslesen

data = np.genfromtxt('data.txt', unpack='True')

data1=ufloat(data[0,:], 0.01)
data2=ufloat(data[0,:], 0.01)

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

def f(x,m, b):
    return m*x+b

#Funktion fitten

guess = [1, 1]
x_plot=np.linspace(0,10, num=1000)
params, covariance = curve_fit(f, unp.nominal_values(data1), unp.nominal_values(data2), p0=guess)

#Daten und Fit plotten

plt.figure()
err1 = unp.std_devs(data1)
err2 = unp.std_devs(data2)

plt.errorbar(unp.nominal_values(data1) + err1, -unp.nominal_values(data2) + err2, fmt='bx', label="Messung")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,10)
plt.ylim(0,10)
plt.xlabel(r'')
plt.ylabel(r'')
plt.savefig('plot1.png')

errors=np.sqrt(np.diag(covariance))

print('m =', params[0], '+/-', errors[0])
print('b =', params[1], '+/-', errors[1])


#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

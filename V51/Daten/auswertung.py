# -*- coding: utf-8 -*-

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

data_linverst = np.genfromtxt('fggverstärker.txt', unpack='True')


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

def g(x,c):
    return c

#Funktion fitten

# guess = [1, 1]
# x_plot=np.linspace(0,10, num=1000)
# params, covariance = curve_fit(f, unp.nominal_values(data1), unp.nominal_values(data2), p0=guess)

#Daten und Fit plotten

amplif = np.ones((19,3))

for i in range(1,4):

    amplif[:,i-1]=data_linverst[2*i-1,:]/data_linverst[2*i,:]
    print(amplif)
    plt.figure()
    err_v = 0.05
    plt.errorbar(data_linverst[0,:], amplif[:,i-1] , fmt='bx', label="Messung zum Widerstand%d" % i)
    #plt.plot(x_plot, f(x_plot, *params), 'r-', label='Fit')
    plt.legend(loc="best", numpoints=1)
    #plt.xlim(0,10)
    #plt.ylim(0,10)
    plt.xlabel(r'Frequenz [kHz]')
    plt.ylabel(r"Verstärkung V' ")
    plt.savefig('../Protokoll/images/plot%d.pdf' % i)

# errors=np.sqrt(np.diag(covariance))
#
# print('m =', params[0], '+/-', errors[0])
# print('b =', params[1], '+/-', errors[1])


#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

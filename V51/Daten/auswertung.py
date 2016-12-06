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

#Daten und Fit plotten

amplif = np.ones((19,4))

for i in range(1,5):

    amplif[:,i-1]=data_linverst[2*i-1,:]/data_linverst[2*i,:]

    x_plot=np.linspace(0.5,10000, num=1000000)
    params_1, covariance_1 = curve_fit(f,np.log(data_linverst[0,13:19]), np.log(amplif[13:19,0]))

    params_2, covariance_2 = curve_fit(f,np.log(data_linverst[0,15:19]), np.log(amplif[15:19,1]))

    params_3, covariance_3 = curve_fit(f,np.log(data_linverst[0,14:19]), np.log(amplif[14:19,2]))

    params_4, covariance_4 = curve_fit(f,np.log(data_linverst[0,10:19]), np.log(amplif[10:19,3]))

    params=np.array([params_1,params_2,params_3,params_4])#+

    print(params)
    plt.figure()
    err_v = 0.05
    if i==4:
        plt.errorbar(np.log(data_linverst[0,0:17]), np.log(amplif[0:17,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)
    else:
        plt.errorbar(np.log(data_linverst[0,:]), np.log(amplif[:,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)

    plt.plot(np.log(x_plot), f(np.log(x_plot), params[:,i-1]), 'r-', label='Fit')
    plt.legend(loc="best", numpoints=1)
    plt.grid()
    plt.ylim(-4,0.1)
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

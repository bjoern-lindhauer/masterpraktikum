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

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

#Daten auslesen

data_linverst = np.genfromtxt('fggverstärker.txt', unpack='True')

data_integ = np.genfromtxt('intfreq.txt', unpack='True')

data_diff = np.genfromtxt('difffreq.txt', unpack='True')

data_schwing_abfall = np.genfromtxt('scope_101_data.txt', unpack='True')

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

def h(x, a, l, k, p):
    return a*np.exp(-x/l)*np.sin(k*x+p)

def h_mod(x,a,k,p):
    return a*np.sin(k*x+p)


#Frequenzen fitten

amplif = np.ones((19,4))

for i in range(1,5):

    amplif[:,i-1]=data_linverst[2*i-1,:]/data_linverst[2*i,:]

    x_plot=np.linspace(0.5,10000, num=1000000)
    params_1, covariance_1 = curve_fit(f,np.log(data_linverst[0,13:19]), np.log(amplif[13:19,0]))

    params_2, covariance_2 = curve_fit(f,np.log(data_linverst[0,15:19]), np.log(amplif[15:19,1]))

    params_3, covariance_3 = curve_fit(f,np.log(data_linverst[0,14:19]), np.log(amplif[14:19,2]))

    params_4, covariance_4 = curve_fit(f,np.log(data_linverst[0,10:17]), np.log(amplif[10:17,3]))

    params=np.array([params_1,params_2,params_3,params_4])

    plt.figure()
    err_v = 0.05
    if i==4:
        plt.errorbar(np.log(data_linverst[0,0:17]), np.log(amplif[0:17,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)
    else:
        plt.errorbar(np.log(data_linverst[0,:]), np.log(amplif[:,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)

    plt.plot(np.log(x_plot), f(np.log(x_plot), params[i-1,0], params[i-1, 1]), 'r-', label='Fit')
    plt.legend(loc="best", numpoints=1)
    plt.grid()
    plt.ylim(-4,0.1)
    plt.xlabel(r'ln(f/kHz)')
    plt.ylabel(r"ln(V')")
    plt.savefig('../Protokoll/images/plot%d.pdf' % i)

#Umkehrintegrator

params_int, covariance_int = curve_fit(f, np.log(data_integ[0,:11]), np.log(data_integ[1,:11]))
x_plot=np.linspace(500,20000, num=100000)

plt.figure()
err_v = 0.05
plt.errorbar(np.log(data_integ[0,:11]), np.log(data_integ[1,:11]) + err_v, fmt='bx', label=r"U$_A$ des Umkehrintegrators")
plt.plot(np.log(x_plot), f(np.log(x_plot), *params_int), 'r-', label='Fit')
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.xlabel(r'ln(f/kHz)')
plt.ylabel(r"ln(U$_A$/V)")
plt.savefig('../Protokoll/images/integrator.pdf')

#Unkehrdifferentiator

params_diff, covariance_diff = curve_fit(f, np.log(data_diff[0,:9]), np.log(data_diff[1,:9]))
x_plot=np.linspace(500,20000, num=100000)

plt.figure()
plt.errorbar(np.log(data_diff[0,:]), np.log(data_diff[1,:]) + err_v, fmt='bx', label=r"U$_A$ des Umkehrdifferentiators")
plt.plot(np.log(x_plot), f(np.log(x_plot), *params_diff), 'r-', label='Fit')
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.xlabel(r'ln(f/kHz)')
plt.ylabel(r"ln(U$_A$/V)")
plt.savefig('../Protokoll/images/differentiator.pdf')

#Abfallende Schwingung

data_schwing_mod = data_schwing_abfall[:,200:1750]
x0=[2.2,0.010,3443,0.038]
x1=[2.5,3500, 0.001]

params_schw, covariance_schw = curve_fit(h, data_schwing_mod[0] , data_schwing_mod[1], p0=x0)
print(params_schw)
x_plot=np.linspace(0,0.018, num=10000)

plt.figure()
plt.plot(data_schwing_mod[0], data_schwing_mod[2], 'bx', label=r'Messwerte')
plt.plot(x_plot, h(x_plot, *params_schw), 'r-', label='nichtlinearer Fit')
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig('../Protokoll/images/schwing_abfall.pdf')

# errors=np.sqrt(np.diag(covariance))
#
# print('m =', params[0], '+/-', errors[0])
# print('b =', params[1], '+/-', errors[1])


#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

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

#T1-Messung

data = np.genfromtxt('T1Messung.txt', unpack='True')

tau=unp.uarray(data[0,:],0)
M=unp.uarray(data[1,:],0.02)

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

def f(t,M0,T1):
    return M0*(1-2*np.exp(-t/T1))

def g(x, b, m):
    return m*x+b

def h(t, M0, D):
    T2=1.227
    G=7.09*10**6
    return M0*np.exp(-2*t/T2)*np.exp(-D*(G**2)*(2*t)**3/12)


#T1-Bestimmung
guess = [1.47, 1]
x_plot=np.linspace(0,10, num=1000)
params, covariance = curve_fit(f, unp.nominal_values(tau), -unp.nominal_values(M), p0=guess)

plt.figure()
errM = unp.std_devs(M)
plt.errorbar(unp.nominal_values(tau), -unp.nominal_values(M) + errM, fmt='bx', label="T1-Messung")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='Nichtlinearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,10)
plt.ylim(0, 2)
plt.xlabel(r'Zeitabstand $\tau$ [s]')
plt.ylabel(r'Magnetisierung M$_z [V]$')
plt.savefig('plotT1.png')

errors=np.sqrt(np.diag(covariance))

print('M0 =', params[0], '+/-', errors[0])
print('T1 =', params[1], '+/-', errors[1])

#T2-Bestimmung

data2 = np.genfromtxt('dataT2a.txt', unpack='True')

tau=unp.uarray(data2[0,:],0)
M=unp.uarray(data2[1,:],0.02)

guess = [1.47, 1]
x_plot=np.linspace(0,1000, num=1000000)
params, covariance = curve_fit(g, unp.nominal_values(tau), np.log(unp.nominal_values(M)), p0=guess)
plt.figure()
errM = unp.std_devs(M)
plt.errorbar(unp.nominal_values(tau), np.log(unp.nominal_values(M)) + errM, fmt='bx', label="T2-Messung")
plt.plot(x_plot, g(x_plot, *params), 'r-', label='Linearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,1010)
plt.ylim(0, 0.8)
plt.xlabel(r'Zeitabstand $\tau$ [$\mu$s]')
plt.ylabel(r'ln(Magnetisierung M$_y$/V)')
plt.savefig('plotT2.png')

errors=np.sqrt(np.diag(covariance))

print('ln(M0) =', params[0], '+/-', errors[0])
print('-(1/T2) =', params[1], '+/-', errors[1])
m0=ufloat(params[0], errors[0])
M0=exp(m0)
print('M0 =', unp.nominal_values(M0) , '+/-', unp.std_devs(M0))
#t1/2-Messung
t= ufloat(params[1],errors[1])
T2=-(1/t)
print('T2=',T2)

t12= ufloat(282,4)
gG=8.8/(4.4*t12)
print(gG)

#Diffusionsskonstante

data2 = np.genfromtxt('DMessung.txt', unpack='True')

tau=unp.uarray(data2[0,:],0)
M=unp.uarray(data2[1,:],0.02)

guess = [500, (2*10**(-7))]
x_plot=np.linspace(0,0.06, num=2000)
params, covariance = curve_fit(h, 2*unp.nominal_values(tau), unp.nominal_values(M), p0=guess)

plt.figure()
errM = unp.std_devs(M)
plt.errorbar(2*unp.nominal_values(tau), unp.nominal_values(M) + errM, fmt='bx', label="D-Messung")
plt.plot(x_plot, h(x_plot, *params), 'r-', label='Nicht-Linearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,0.05)
plt.ylim(0, 570)
plt.xlabel(r'Zeitabstand 2$\tau$ [$\mu$s]')
plt.ylabel(r'Magnetisierung M$_y$ [V]')
plt.savefig('plotD.png')

errors=np.sqrt(np.diag(covariance))

print('M0 =', params[0], '+/-', errors[0])
print('D =', params[1], '+/-', errors[1])

#Viskosität
t=947
de=0.45
rho=1000
a=1.024*10**(-9)
eta=rho*a*(t-de)
print(eta)

#Molekülradius
k=1.38*10**(-23)
r = (k*293)/(6*np.pi*1.11*10**(-7)*0.97*10**(-3))
print(r)
#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

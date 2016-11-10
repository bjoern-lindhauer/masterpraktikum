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
from uncertainties.unumpy import log10,log,exp,sqrt,sin,arctan
import sympy

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

#Daten auslesen

#Daten zur Bestimmung der Wellenlänge formartieren
wave = np.genfromtxt('wavelength.txt', unpack='True')
n=wave[0,:]
d=unp.uarray(wave[1,:]*10**(-2),0.3)

#Polarisationsmessdaten einlesen
polarisation = np.genfromtxt('polarisation.txt', unpack='True')

#Daten zur Vermessung der TEM00-Mode einlesen
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
    return I0*np.exp(-2*(r-r0)**2/w**2)

def Lambda(d,L,n,g):
    return (sin(arctan(d/L)))/(abs(n)*g)

def Ip(phi,I0,delta):
    return I0*np.cos(phi+delta)**2

def f(x,r):
    return (1-x/r)

def g(x,r1,r2):
    return (1-x/r1)*(1-x/r2)

#Konstanten defininieren

g1=100*10**(3)
L=ufloat(140.5*10**(-2),0.1)

#Wellenlängebestimmung

l=Lambda(d,L,n,g1)
l_mean=l.mean();
print(l_mean)


#Polarisation fitten

params_pol, covariance_pol = curve_fit(Ip, np.deg2rad(polarisation[0]), polarisation[1])
x_plot=np.linspace(0,180,num=1000)

plt.figure()
plt.plot(x_plot, Ip(np.deg2rad(x_plot), *params_pol), 'r-', label='Nichtlinearer Fit')
plt.plot(polarisation[0], polarisation[1], 'gx', label='Messdaten Polarisation')
plt.xlabel(r'Winkel $\phi$')
plt.ylabel(r'Intensität/$\mu$A')
plt.legend(loc="best", numpoints=1)
plt.savefig('../Protokoll/images/polarisaton.pdf')
plt.close()

#Vorbereitungsaufgabe

x_plot=np.linspace(0,10, num=1000)
plt.figure()
plt.plot(x_plot, f(x_plot, 1.4), 'r-', label='Planarer und r=1m Spiegel')
plt.plot(x_plot, g(x_plot, 1, 1.4), 'gx', label='r=1m Spiegel und r=1.4m Spiegel')
plt.legend(loc="best", numpoints=1)
plt.xlim(0,2)
plt.ylim(0,1)
plt.xlabel(r'L/m')
plt.ylabel(r'g$_1*$g$_2$')
plt.savefig('../Protokoll/images/vorbereitung.pdf')
plt.close()

#TEM00-Mode fitten

x0=[0.1,2,4]
params_tem, covariance_tem = curve_fit(I, tem00[0], tem00[1], p0=x0)

x_plot=np.linspace(-15,15,num=1000)


plt.figure()
plt.plot(x_plot, I(x_plot, *params_tem), 'r-', label='Nichtlinearer Fit')
plt.plot(tem00[0], tem00[1], 'gx', label='Messdaten TEM00-Mode')
plt.xlabel('r/mm')
plt.ylabel('I/$\mu$A')
plt.legend(loc="best", numpoints=1)
plt.savefig('../Protokoll/images/tem00.pdf')
plt.close()

# errors=np.sqrt(np.diag(covariance))
#
# print('m =', params[0], '+/-', errors[0])
# print('b =', params[1], '+/-', errors[1])
#

#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

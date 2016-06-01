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

data = np.genfromtxt('gleichrichter.txt', unpack='True')

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

def f(t,a,b):
    return a*np.cos(t+b)
    
def g(t):
    return t*9.58*10**(-3)*2*np.pi
    
def h(x,y):
    return 2*10**((x-y)/20)
    
def kontrast(x2, x1):
    return (x2-x1)/(x2+x1)
    
def j(dt, ft):
    return ((-1/(dt*ft))+unp.sqrt((1/(dt*ft))**2+1))


#Fit für Gleichrichter
x_plot=np.linspace(0,7, num=1000)
params, covariance = curve_fit(f, g(data[0,:]) , data[1,:])

plt.figure()
plt.plot(g(data[0,:]), data[1,:], 'bx', label="Gleichrichter-Messung")
plt.plot(x_plot, f(x_plot, *params), 'r-', label='Nichtlinearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlim(0, 7)
plt.ylim(-0.5, 0.5)
plt.xlabel(r'Phasendifferenz $\Phi$')
plt.ylabel(r'Signalstaerke U [V]')
plt.savefig('plotgleich.png')

errors=np.sqrt(np.diag(covariance))

print('a =', params[0], '+/-', errors[0])
print('b =', params[1], '+/-', errors[1])


#Modulationsgrad berechnen
Umax=ufloat(0.290,0.05)
Umin=ufloat(0.160,0.05)

print(kontrast(Umax,Umin))

Ut=ufloat(-26.12, 0.01)
Um=ufloat(-48.67, 0.01)
print('m=', h(Um,Ut))

dt=ufloat(375*10**(-9),25*10**(-9))

print('m=', j(dt, 10**6))

#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

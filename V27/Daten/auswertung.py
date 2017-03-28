# -*- coding: utf-8 -*-

import numpy as np
from numpy.lib import scimath
import scipy
from scipy.optimize import curve_fit
from scipy import constants
from scipy.stats import sem
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import log10,log,exp,sqrt
import sympy
import cmath

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

h=constants.value(u'Planck constant')
c=constants.value(u'speed of light in vacuum')
mb = constants.value(u'Bohr magneton')
l_red = 480*10**(-9)
l_blue = 643.8*10**(-9)
ld_red = 48.9*10**(-12)
ld_blue = 26.9*10**(-12)

data_b = np.genfromtxt('magnetfeld.txt', unpack='True')

data_deltas_red = np.genfromtxt('deltas_red.txt', unpack='True', skip_header=1)
data_deltas_blue = np.genfromtxt('deltas_blue.txt', unpack='True', skip_header=1)

def f(x,a,b):
    return a*x+b

def dlambda(ds, Ds, Dl):
    return 0.5*ds/Ds*Dl

def lande(l,B,dl):
    return h*c/(l**2*mb*B)*dl

params, covariance = curve_fit(f,data_b[0],data_b[1])
errors = np.sqrt(np.diag(covariance))

params_err = unp.uarray(params,errors)

print('m= {:.2u}'.format(params_err[0]))
print('b= {:.2u}'.format(params_err[1]))
x_plot = np.linspace(0,20, 1000)

plt.figure()
plt.plot(data_b[0], data_b[1], 'bx', label="gemessenes Magnetfeld")
plt.plot(x_plot, f(x_plot, *params), 'r-', label="linearer Fit")
plt.xlabel(r"Stromstärke I / A")
plt.ylabel(r"Magnetfeld B / T")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("../Protokoll/images/magnetfeld.pdf")

#calculation of change of wave-length

delta_l_red = dlambda(data_deltas_red[1],data_deltas_red[0],ld_red)

delta_l_blue = np.ones((2,10))

delta_l_blue[0] = dlambda(data_deltas_blue[1], data_deltas_blue[0], ld_blue)
delta_l_blue[1] = dlambda(data_deltas_blue[2], data_deltas_blue[0], ld_blue)

dl_red = ufloat(np.mean(delta_l_red), sem(delta_l_red))
dl_blue = unp.uarray([np.mean(delta_l_blue[0]),np.mean(delta_l_blue[1])] , [sem(delta_l_blue[0]),sem(delta_l_blue[1])])

print('dl(643.8nm) = {:.2u}'.format(dl_red))
print('B = {:.5}'.format(f(11.5,*params)))
print('dl(480.0nm) = {:.2u}'.format(dl_blue[0]))
print('B = {:.5}'.format(f(5.5,*params)))
print('dl(480.0nm) = {:.2u}'.format(dl_blue[1]))
print('B = {:.5}'.format(f(18,*params)))

# for i in range (0,len(delta_l_red)):
#     print(i, "\t", '{:.4}'.format(delta_l_red[i]), "\n")
#
# for i in range (0,len(delta_l_blue[0])):
#     print(i, "\t", '{:.4}'.format(delta_l_blue[0,i]), "\n")
#
# for i in range (0,len(delta_l_blue[1])):
#     print(i, "\t", '{:.4}'.format(delta_l_blue[1,i]), "\n")

print('Die Landé-Faktoren lauten:')
print('\t', 'rot:', lande(l_red,f(11.5,*params),dl_red))
print('\t', 'blau, sigma:', lande(l_blue,f(5.5,*params),dl_blue[0]))
print('\t', 'blau, pi:', lande(l_blue,f(18,*params),dl_blue[1]))

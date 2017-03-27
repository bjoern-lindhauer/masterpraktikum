# -*- coding: utf-8 -*-

import numpy as np
from numpy.lib import scimath
import scipy
from scipy.optimize import curve_fit
from scipy import constants
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import log10,log,exp,sqrt
import sympy
import cmath

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

data = np.genfromtxt('magnetfeld.txt', unpack='True')

def f(x,a,b):
    return a*x+b

params, covariance = curve_fit(f,data[0],data[1])
errors = np.sqrt(np.diag(covariance))

params_err = unp.uarray(params,errors)

print('m= {:.2u}'.format(params_err[0]))
print('b= {:.2u}'.format(params_err[1]))
x_plot = np.linspace(0,20, 1000)

plt.figure()
plt.plot(data[0], data[1], 'bx', label="gemessenes Magnetfeld")
plt.plot(x_plot, f(x_plot, *params), 'r-', label="linearer Fit")
plt.xlabel(r"Stromstärke I / A")
plt.ylabel(r"Magnetfeld B / T")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("../Protokoll/images/magnetfeld.pdf")

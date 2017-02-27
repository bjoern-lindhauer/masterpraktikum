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

l=1.54*10**(-10) #wavelength in Angstrom
k=2*np.pi/l #wavenumber

ai = np.linspace(0,5,10000)*np.pi/180 #incident angle
qz = 4*np.pi/(l*10**(10))*np.sin(ai) #change of wave vector


def reflectrometry(a,n1,n2,n3,s1,s2):

    #z-components
    kz1=k*np.sqrt(n1^2-np.square(np.cos(a)))
    kz2=k*np.sqrt(n2^2-np.square(np.cos(a)))
    kz3=k*np.sqrt(n3^2-np.square(np.cos(a)))

    #modified Fresnel-coefficients

    r12=(kz1-kz2)/(kz1+kz2)*np.exp(-2*kz1*kz2*s1^2)
    r23=(kz2-kz3)/(kz2+kz3)*np.exp(-2*kz2*kz3*s2^2)
    x2=exp(-2*complex(0,1)*kz2*z2)*r23
    x1=(r12+x2)/(1+r12*x2)

    return x1




plt.figure()

plt.yscale("log")
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.savefig("reflectrometry.pdf")

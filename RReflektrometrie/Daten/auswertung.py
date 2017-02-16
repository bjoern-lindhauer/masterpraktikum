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

l=1.54 #wavelength in Angstrom
k=2*np.pi/l #wavenumber

def kz(n,a):
    return(k*(n**2-np.cos(a)**2)**(1/2))

def  rji(k,nj,ni,aj,ai):
    kj=kz(nj,aj)
    ki=kz(nj,aj)
    return((kj-ki)/(kj+ki))

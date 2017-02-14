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

data_phase = np.genfromtxt('freqphasggverst.txt', unpack='True')

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

#Defining functions

def f(x,m,b):
    return m*x+b

def g(x,c):
    return c

def h(x,a,l):
    return a*np.exp(-x/l)

def const(x,c):
    return c*np.ones(np.shape(x))


#fitting frequencies

amplif = np.ones((19,4))
vgrenz= np.ones(5)
vgrenz_errors= unp.uarray(np.ones(4),0.01)
th_l = np.array([15,15,15,11]) #lower threshold for fit
th_u = np.array([20,20,20,18]) #higher threshold for fit
params = [np.zeros((2,1)),np.zeros((2,1)),np.zeros((2,1)),np.zeros((2,1))]
covariance = [np.zeros((2,2)),np.zeros((2,2)),np.zeros((2,1)),np.zeros((2,2))]
errors_cov = [np.zeros((2,1)),np.zeros((2,1)),np.zeros((2,1)),np.zeros((2,1))]

for i in range(1,5):

    amplif[:,i-1]=data_linverst[2*i-1,:]/data_linverst[2*i,:]

    x_plot=np.linspace(0.5,10000, num=1000000)
    params[i-1], covariance[i-1] = curve_fit(f,np.log(data_linverst[0,th_l[i-1]-1:th_u[i-1]]),np.log(amplif[th_l[i-1]-1:th_u[i-1],i-1]))
    errors_cov[i-1]=np.sqrt(np.diag(covariance[i-1]))

    plt.figure()
    err_v = 0.05
    plt.errorbar(data_linverst[0,0:th_l[i-1]-1], amplif[0:th_l[i-1]-1,i-1] + err_v , fmt='bx', label="Nicht für den Fit verwendet")
    plt.errorbar(data_linverst[0,th_l[i-1]:th_u[i-1]], amplif[th_l[i-1]:th_u[i-1],i-1] + err_v , fmt='yx', label="Für den Fit verwendet")

    vgrenz[i-1] = np.mean(np.log(amplif[0:3,i-1])-0.5*np.log(2))
    vgrenz_errors[i-1] = ufloat(vgrenz[i-1],np.std(np.log(amplif[0:3,i-1])-0.5*np.log(2)))

    plt.plot(x_plot, np.exp(f(np.log(x_plot), *params[i-1])), 'r-', label='Fit')
    plt.plot([5*10**(-2),5*10**4], const([5*10**(-2),5*10**4], np.exp(vgrenz[i-1])), 'g-', label=r"Grenzwert zur Grenzfrequenz")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="lower left", numpoints=1)
    plt.grid()
    plt.xlim(5*10**(-2),10**5)
    plt.ylim(10**(-2),2)
    plt.xlabel(r'f/kHz')
    plt.ylabel(r"V'")
    name = "../Protokoll/images/plot%d_varied.pdf" % i
    plt.savefig(name)
    plt.close()


print('Die Steigungen lauten:')
print('m_1 =', '%.3f' % params[0][0], '+/-', '%.3f' % errors_cov[0][0])
print('m_2 =', '%.3f' % params[1][0], '+/-', '%.3f' % errors_cov[1][0])
print('m_3 =', '%.3f' % params[2][0], '+/-', '%.3f' % errors_cov[2][0])
print('m_4 =', '%.3f' % params[3][0], '+/-', '%.3f' % errors_cov[3][0])

m_1,m_2,m_3,m_4 =ufloat(params[0][0], errors_cov[0][0]), ufloat(params[1][0], errors_cov[1][0]), ufloat(params[2][0], errors_cov[2][0]), ufloat(params[3][0], errors_cov[3][0])

b_1,b_2,b_3,b_4 =ufloat(params[0][1], errors_cov[0][1]), ufloat(params[1][1], errors_cov[1][1]), ufloat(params[2][1], errors_cov[2][1]), ufloat(params[3][1], errors_cov[3][1])

print('Die Grenzfrequenzen lauten:')
fgrenz = np.array([exp((vgrenz[0]-b_1)/m_1),exp((vgrenz[1]-b_2)/m_2),exp((vgrenz[2]-b_3)/m_3),exp((vgrenz[3]-b_4)/m_4)])
print(fgrenz)

print('Die Verstärkungen lauten:')
vgrenz_exp = exp(vgrenz_errors)*sqrt(2)
print(vgrenz_exp)

vtheor = np.array([1.0,0.1,0.57,100])

print('Die Leerlaufverstärkungen lauten:')
vleer = (vgrenz_exp*vtheor)/(vtheor-vgrenz_exp)
print(vleer)

print('Das Verstärkung-Bandbreite-Produkt lautet:')
bbprodukt = vgrenz_exp*fgrenz
print(bbprodukt)


#Umkehrintegrator

params_int, covariance_int = curve_fit(f, np.log(data_integ[0,:11]), np.log(data_integ[1,:11]))
errors_int = np.sqrt(np.diag(covariance_int))

print(params_int, errors_int)

x_plot=np.linspace(500,20000, num=100000)

plt.figure()
plt.plot(data_integ[0,:], data_integ[1,:], 'bx', label=r"U$_A$ des Umkehrintegrators")
plt.plot(x_plot, np.exp(f(np.log(x_plot), *params_int)), 'r-', label='Fit')
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.xlabel(r'f/kHz')
plt.ylabel(r"U$_A$/V")
plt.xscale("log")
plt.yscale("log")
plt.xlim(40,2*10**4)
plt.ylim(10**(-1),4)
plt.savefig('../Protokoll/images/integrator.pdf')
#Unkehrdifferentiator

params_diff, covariance_diff = curve_fit(f, np.log(data_diff[0,:9]), np.log(data_diff[1,:9]))
errors_diff = np.sqrt(np.diag(covariance_diff))

print(params_diff, errors_diff)

x_plot=np.linspace(500,20000, num=100000)

plt.figure()
plt.errorbar(data_diff[0,:], data_diff[1,:] + err_v, fmt='bx', label=r"U$_A$ des Umkehrdifferentiators")
plt.plot(x_plot, np.exp(f(np.log(x_plot), *params_diff)), 'r-', label='Fit')
plt.legend(loc="best", numpoints=1)
plt.grid()
plt.xlabel(r'f/kHz')
plt.ylabel(r"U$_A$/V")
plt.xscale("log")
plt.yscale("log")
plt.xlim(9,2*10**4)
plt.ylim(10**(-1),4)
plt.savefig('../Protokoll/images/differentiator.pdf')

#Abfallende Schwingung

data_schwing_mod = data_schwing_abfall[:,200:1750]

#Peak-Detection

peaks = np.ones((2,16))

for i in range (1,16):

     j = (i-1)*100
     k = j+100
     index = np.argmax(data_schwing_mod[2,j:k])
     index=index+j
     peaks[0,i] = data_schwing_mod[0,index]
     peaks[1,i] = data_schwing_mod[2,index]

peaks = np.delete(peaks, 0, 1)

print(peaks)
x0=[3,0.008]

params_exp, covariance_exp = curve_fit(h, peaks[0], peaks[1], p0=x0)
errors_exp = np.sqrt(np.diag(covariance_exp))
T = ufloat(params_exp[1], errors_exp[1])
print('Die Abklindauer beträgt:')
print(T)


x_plot=np.linspace(0,0.018, num=10000)

plt.figure()
plt.plot(data_schwing_mod[0], data_schwing_mod[2], 'bx', label=r'Messwerte')
plt.plot(x_plot, h(x_plot, *params_exp), 'r-', label='nichtlinearer Fit')
plt.legend(loc="best", numpoints=1)
plt.xlabel('t/s')
plt.ylabel(r'Spannung $U_A$/V')
plt.grid()
plt.savefig('../Protokoll/images/schwing_abfall.pdf')

#Phasebeziehung Gegenverstärker

plt.figure()
plt.plot(data_phase[0], np.deg2rad(data_phase[1]), 'bx', label=r'Messwerte')
plt.yticks([np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
          [r'$+\pi/2$', r'$+\pi$', r'$+3\pi/2$', r'$2\pi$'])
plt.xlabel(r'Frequenz f/Hz')
plt.ylabel(r'Phase $\phi$')
plt.legend(loc='best', numpoints=1)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.grid()
plt.savefig('../Protokoll/images/phase_frequenz.pdf')


#x1,x2, x3 = sympy.var('M_{z} M_{0} \tau')
#f = -x3/(sympy.log((x1-x2)/(2*x2)))
#print(error(f))

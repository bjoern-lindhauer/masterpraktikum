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

#Funktionen definieren

def f(x,m,b):
    return m*x+b

def g(x,c):
    return c

def h(x,a,l):
    return a*np.exp(-x/l)

def const(x,c):
    return c*np.ones(np.shape(x))


#Frequenzen fitten

amplif = np.ones((19,4))
vgrenz= np.ones(5)
vgrenz_errors= unp.uarray(np.ones(4),0.01)

for i in range(1,5):

    amplif[:,i-1]=data_linverst[2*i-1,:]/data_linverst[2*i,:]

    x_plot=np.linspace(0.5,10000, num=1000000)
    params_1, covariance_1 = curve_fit(f,np.log(data_linverst[0,13:19]), np.log(amplif[13:19,0]))
    errors_1=np.sqrt(np.diag(covariance_1))

    params_2, covariance_2 = curve_fit(f,np.log(data_linverst[0,15:19]), np.log(amplif[15:19,1]))
    errors_2=np.sqrt(np.diag(covariance_2))

    params_3, covariance_3 = curve_fit(f,np.log(data_linverst[0,14:19]), np.log(amplif[14:19,2]))
    errors_3=np.sqrt(np.diag(covariance_3))

    params_4, covariance_4 = curve_fit(f,np.log(data_linverst[0,10:17]), np.log(amplif[10:17,3]))
    errors_4=np.sqrt(np.diag(covariance_4))

    params=np.array([params_1,params_2,params_3,params_4])

    plt.figure()
    err_v = 0.05
    if i==4:
        plt.errorbar(np.log(data_linverst[0,0:17]), np.log(amplif[0:17,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)
    else:
        plt.errorbar(np.log(data_linverst[0,:]), np.log(amplif[:,i-1]) + err_v , fmt='bx', label="Messung zum Widerstand %d" % i)

    vgrenz[i-1] = np.mean(np.log(amplif[0:3,i-1])-0.5*np.log(2))
    vgrenz_errors[i-1] = ufloat(vgrenz[i-1],np.std(np.log(amplif[0:3,i-1])-0.5*np.log(2)))

    plt.plot(np.log(x_plot), f(np.log(x_plot), params[i-1,0], params[i-1, 1]), 'r-', label='Fit')
    plt.plot([-4,10], const([-4,10], vgrenz[i-1]), 'g-', label=r"Grenzwert zur Grenzfrequenz")
    plt.legend(loc="best", numpoints=1)
    plt.grid()
    plt.ylim(-4,0.1)
    plt.xlabel(r'ln(f/kHz)')
    plt.ylabel(r"ln(V')")
    name = "../Protokoll/images/plot%d.pdf" % i
    plt.savefig(name)


print('Die Steigungen lauten:')
print('m_1 =', '%.3f' % params_1[0], '+/-', '%.3f' % errors_1[0])
print('m_2 =', '%.3f' % params_2[0], '+/-', '%.3f' % errors_2[0])
print('m_3 =', '%.3f' % params_3[0], '+/-', '%.3f' % errors_3[0])
print('m_4 =', '%.3f' % params_4[0], '+/-', '%.3f' % errors_4[0])

m_1 =ufloat(params_1[0], errors_1[0])
m_2 =ufloat(params_2[0], errors_2[0])
m_3 =ufloat(params_3[0], errors_3[0])
m_4 =ufloat(params_4[0], errors_4[0])

b_1 =ufloat(params_1[1], errors_1[1])
b_2 =ufloat(params_2[1], errors_2[1])
b_3 =ufloat(params_3[1], errors_3[1])
b_4 =ufloat(params_4[1], errors_4[1])

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
errors_diff = np.sqrt(np.diag(covariance_diff))

print(params_diff, errors_diff)

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

x0=[3,0.008]

params_exp, covariance_exp = curve_fit(h, peaks[0], peaks[1], p0=x0)
x_plot=np.linspace(0,0.018, num=10000)

plt.figure()
plt.plot(data_schwing_mod[0], data_schwing_mod[2], 'bx', label=r'Messwerte')
plt.plot(x_plot, h(x_plot, *params_exp), 'r-', label='nichtlinearer Fit')
plt.legend(loc="best", numpoints=1)
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

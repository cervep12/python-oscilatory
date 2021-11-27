# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 17:58:50 2021

@author: Petr
"""
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Model
from min_ctverce import *
from scipy.misc import derivative
from zfmvyp import *
from scipy import stats
from chyba_aritm_prum import *
import pandas as pd
import lmfit
from lmfit import Model
def fit_function(t,A,B,C,D,E):
    return C * np.exp(-D*t)*np.sin(A*t + B) + E
def rezon(omeg,omeg0,delta,B):
    return B/((omeg0**2-omeg**2)**2 + 4*delta**2*omeg**2)**0.5
def rezon_netl(omeg,omeg0,B):
    return B/np.abs(omeg0**2-omeg**2)



kalibrU = [8.04,4,8,12,16,20]
omU = np.array([32.4,13.8,29.3,44.3,60.4,76.5])/60
T = np.array([1.68,1.67,1.67,1.67,1.67,1.68,1.67,1.67,1.67,1.67])
chybaT = np.array([0.001,0.0013,0.0014,0.0011,0.0012,0.0011,0.0027,0.0038,0.0025,0.0031,])
print("chybaT, prumT: ", chyba_aritm_prum(T),np.mean(T))
v = komb_mereni(chybaT,T)
print("zahrnu nepres: ", v)
omega0 = 2*np.pi/v[0]
sigm_om0 = omega0/v[0]*v[1]
print("om0 , chyba om0: ", omega0, sigm_om0)
print("------------------------")
#print(np.mean(omega0),np.mean(omega0)/)



Atl = np.array([0.284,0.286,0.793,1.076,1.505,0.611,0.377,0.24])
Anetl = np.array([0.193,0.38,0.763,2.16,1.203,0.735,0.339,0.286])
Utl =  np.array([7.7,8.15,9.05,9.4,9.75,10.4,11,11.75])*2*math.pi*0.0631
Unetl =  np.array([7.3,8.2,8.95,9.4,10.05,10.4,11.02,11.75])*2*math.pi*0.0631
ftl = np.array([28.4,30.2,33.3,34.5,35.8,38.7,41.0,44.1])/60*2*math.pi
fnetl = np.array([26.9,30.3,33.2,34.2,37.3,38.8,41.4,44.1])/60*2*math.pi

plt.figure(figsize = (6.5,4),dpi=300)
ptl = [3.75,0.15,1.6]
plt.scatter(ftl,Atl,color = "purple",marker = "x" ,label = "tlumené kmity")
popttl,pcovtl = curve_fit(rezon,ftl,Atl,p0=ptl)
sigtl = np.sqrt(np.diag(pcovtl))
a = np.linspace(2.6,5,1001)
fit_hodnoty = rezon(a,*popttl)
plt.plot(a,fit_hodnoty, label = r"$\varphi_m$ = $h_{tl}(\Omega)$")
#--------------------------------------------


pnetl = [3.7,0.1,1.759]
plt.scatter(fnetl,Anetl,color = "black",marker = "D",s = 20, label = "netlumené kmity")
poptnetl,pcovnetl = curve_fit(rezon,fnetl,Anetl,p0=pnetl)
signetl = np.sqrt(np.diag(pcovnetl))
a = np.linspace(2.6,3.67,1000)
b = np.linspace(3.73,5,1000)
fit_hodnoty2 = rezon(a,*poptnetl)
fit_hodnoty3 = rezon(b,*poptnetl)
plt.plot(a,fit_hodnoty2,color = "red",ls = "--", label = r"$\varphi_m$ = $h_{netl}(\Omega)$")
plt.plot(b,fit_hodnoty3,color = "red",ls = "--")
plt.axvline(x = 3.700,color = "red",ls = ":",label = "$\Omega = 3.701$")
plt.axvline(x = 3.755,color = "blue",ls = ":",label = "$\Omega = 3.755$")
plt.grid()
plt.legend(fontsize = 12)
plt.xlabel(r"$\Omega[rad\cdot s^{-1}]$", fontsize = 13)
plt.ylabel(r"$\varphi_m[rad]$", fontsize = 13)
plt.savefig("rezon2.eps")
plt.show()
print("netl: ", poptnetl,signetl)
print("tl: ", popttl,sigtl)
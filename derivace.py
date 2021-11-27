# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 14:50:12 2021

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
data = pd.read_csv("C:\\Users\\Petr\\Desktop\\pohl\\p1.5.txt",delimiter = "\t",names = ['A','B'])
cas1= data["A"].to_numpy()
uhel1 =  data["B"].to_numpy()-0.87
data = pd.read_csv("C:\\Users\\Petr\\Desktop\\pohl\\p2.4.txt",delimiter = "\t",names = ['A','B'])
cas2= data["A"].to_numpy()
uhel2 =  data["B"].to_numpy()
p03 = [1,1,1,1,1]
popt3,pcov3 = curve_fit(fit_function,cas1,uhel1,p0=p03)


func1 = lambda t:popt3[2] * np.exp(-popt3[3]*t)*np.sin(popt3[0]*t + popt3[1]) + popt3[4]

sig3 = np.sqrt(np.diag(pcov3))
print(popt3)
print("--------")
def deriv(t):
    return derivative(func1,t)
der = deriv(cas1)
y = []
for i in range(len(cas1)):
    y.append(1.353)

plt.figure(figsize = (6.5,4),dpi=150)

x1 = plt.plot(cas1,uhel1,color = "black",ls = "--",label = "netlumené kmity")
x2 = plt.plot(cas2,uhel2,color = "red",label = "tlumené kmity")
#plt.plot(uhel,der)
# plt.plot(cas1,uhel,label = "naměřené hodnoty",color = "black")
# plt.plot(cas1,np.array(y),label = r"$\varphi = 1.353 \cdot t$",color = "red")
plt.xlabel(r"$t{[s]}$",size = 14)
plt.ylabel(r"$\varphi[rad]$",size = 14)
# plt.legend(fontsize = 12)
plt.grid(color = "whitesmoke")
plt.legend(fontsize = 13)
plt.savefig("kmity.eps")
plt.show()
#plt.plot(uhel,der)
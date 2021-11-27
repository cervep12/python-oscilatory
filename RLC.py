# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 10:05:59 2021

@author: Petr
"""

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
def rezon(f,f0,Im,Q):
    return Im/((f0/f-f/f0)**2*Q**2+ 1)**0.5

fpr = [217.4,211.9,204.9,208.3,208.3,210.1,208.3,208.3,206.6,210.1]
Cx = [403,381,373,390,390]
Csig = [1,1,1,1,1]



data = pd.read_csv("C:\\Users\\Petr\\Desktop\\rezon.txt",delimiter = "\t",names = ['A','B',"C"])
f= np.array(data["A"].to_numpy())
I = np.array(data["B"].to_numpy())/5*10**3
Ij = np.array(data["C"].to_numpy())/5*10**3




plt.figure(figsize = (6.8,4.3),dpi=300)
ptl = [205,35,3]
plt.scatter(f,I,color = "black",s=20 ,label = "cívka " + "\u0332".join("bez ") +"jádra")
popttl,pcovtl = curve_fit(rezon,f,I,p0=ptl)
sigtl = np.sqrt(np.diag(pcovtl))
a = np.linspace(135,275,1001)
fit_hodnoty = rezon(a,*popttl)
plt.plot(a,fit_hodnoty, label = r"$I_0$ = $h_{bez}(f)$")
#--------------------------------------------


pnetl = [205,22,1000]
plt.scatter(f,Ij,color = "purple",marker = "x",s = 40, label = "cívka "+ "s"+"\u0332"+" jádrem")
poptnetl,pcovnetl = curve_fit(rezon,f,Ij,p0=pnetl)
signetl = np.sqrt(np.diag(pcovnetl))
a = np.linspace(135,275,1000)
fit_hodnoty2 = rezon(a,*poptnetl)
plt.plot(a,fit_hodnoty2,color = "red",ls = "--", label = r"$I_0$ = $h_{s}(f)$")

plt.axvline(x = 208.159,color = "red",ls = ":",label = "$f = 208$")
plt.axvline(x = 207,color = "blue",ls = ":",label = "$f = 207$")
plt.grid()
plt.legend(fontsize = 12)
plt.xlabel(r"$f~[kHz]$", fontsize = 13)
plt.ylabel(r"$I_0[mA]$", fontsize = 13)
plt.savefig("RLC.eps")
plt.show()
print("netl: ", poptnetl,signetl)
print("tl: ", popttl,sigtl)
print(np.mean(fpr),chyba_aritm_prum(fpr))

print(np.mean(Cx),chyba_aritm_prum(Cx))
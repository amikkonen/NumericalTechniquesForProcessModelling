# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 23:48:11 2017

@author: ami
"""
import numpy as np

rho = 1000
L   = 10
d   = 0.05
nu  = 1e-6


Vs=[1,2,3,4,5,6,7,8,9,10]

dps = []
for V in Vs:
    Re = V*d/nu
    f = (-1.8*np.log10(6.9/Re))**-2
    dp = 0.5*rho*V**2*L/d*f
    
    dps.append(dp)
    
    print("V  =", V)    
    print("Re =", Re)
    print("f  =", f)
    print("dp =", dp, "Pa")
    print()

print(dps)






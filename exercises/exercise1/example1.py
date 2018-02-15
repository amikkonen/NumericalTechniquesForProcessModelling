# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 23:48:11 2017

@author: ami
"""
import scipy as sp
from matplotlib import pyplot as plt

def solver(rho, L, d, nu, V):
    Re = V*d/nu
    f  = (-1.8*sp.log(6.9/Re))**-2
    dp = 0.5*rho*V**2*L/d*f
    
    print("V  =", V)    
    print("Re =", Re)
    print("f  =", f)
    print("dp =", dp, "Pa")
    
    return dp

rho = 1000
L   = 10
d   = 0.05
nu  = 1e-6
V   = sp.array([1,2,3,4,5,6,7,8,9,10])

dp = solver(rho, L, d, nu, V)

plt.plot(V, dp)
plt.xlabel("V (m/s)")
plt.ylabel("dp (Pa)")
plt.grid()


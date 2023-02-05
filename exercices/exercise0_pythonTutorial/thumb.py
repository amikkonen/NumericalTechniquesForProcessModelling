#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:32:18 2021

@author: ami
"""

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI

def solver(rho, L, d, nu, V):
    Re = V*d/nu
    f  = (-1.8*np.log10(6.9/Re))**-2
    dp = 0.5*rho*V**2*L/d*f
    return dp    
    
if __name__ == "__main__":    
        
    T = 1000
    rho = PropsSI("D", "T", T, "P", 1e6, "air")
    L   = 10
    d   = 0.05
    d2  = 0.025
    mu  = PropsSI("V", "T", T, "P", 1e6, "air")
    nu = mu / rho
    
    V = np.array([1,2,3,4,5,6,7,8,9,10])
    
    dp = solver(rho, L, d, nu, V)
    dp2 = solver(rho, L, d2, nu, V)
    
    plt.plot(V, dp)
    plt.plot(V, dp2)
    plt.xlabel('V (m/s)')
    plt.ylabel('dp (Pa)')
    plt.grid()
    plt.xlim(1,10)


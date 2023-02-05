# -*- coding: utf-8 -*-
"""
@author: Antti Mikkonen, a.mikkonen@iki.fi, 2019
"""

import scipy as sp
from matplotlib import pyplot as plt
from scipy.integrate import solve_bvp


C2K = 273.15
# Stefan-Boltman constant
sigma = 5.67036713e-8

def analytical(x,L,T_B, T_inf, h, kf, d):
    
    beta = (4*h/kf/d)**0.5
    return (T_B-T_inf)*(sp.cosh(beta*(L-x))/sp.cosh(beta*L) ) + T_inf


def numerical(x_ar,T_B, T_inf, h, kf, d, emissivity):
        
#    def ode(x, y):
#        return sp.array([, # place holder
#                          # place holder
#                         ])
#
#    def bc(ya, yb):
#        return sp.array([, # place holder
#                         # place holder
#                        ]) 
#    
#    y_init = [sp.zeros_like(x_ar), # place holder
#              sp.zeros_like(x_ar) # place holder
#              ]
#
#    sol = solve_bvp(ode, bc, x_ar, y_init)
##    print(sol)
#    assert sol.success

    return sp.zeros_like(x_ar)+T_B # place holder


#############################################################################
T_B   = 60+C2K
T_inf = 20+C2K
h     = 5
d     = 0.001
L     = 0.1
kf    = 237
emissivity = 0.0

# Numerical solution locations
n = 10
x_numerical = sp.linspace(0, L, n)

# Analytical solution
x_analytical = sp.linspace(0, L, 1000)
T_analytical = analytical(x_analytical,L,T_B, T_inf, h, kf, d)

# Numerical solution
T_numerical = numerical(x_numerical,T_B, T_inf, h, kf, d, emissivity)

# Plotting
plt.plot(x_analytical/L, T_analytical-C2K, label="analytical")
plt.plot(x_numerical/L, T_numerical-C2K, "d", label="numerical")

plt.xlim(0,1)
plt.legend(frameon=False)
plt.xlabel("x/L")
plt.ylabel("T (C)")
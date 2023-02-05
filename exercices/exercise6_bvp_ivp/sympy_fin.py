#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 11:18:58 2019

@author: ami
"""

import sympy as sm
import numpy as np
from matplotlib import pyplot as plt



# Define analytical symbols
h, k, d, x, T_inf, T_B = sm.symbols("h k d x T_inf T_")
T = sm.Function('T')

# Define ode and intial value
#ode = sm.Derivative(T(x),x,x) + a*(T(x)-T_inf)
ode = sm.Derivative(T(x),x,x) + 4*h/k/d*(T(x)-T_inf) 
print(ode)


intode = sm.dsolve(sm.integrate(ode, x))
print(intode)




#ode1 = sm.dsolve(ode)
#print(ode1)



#def get_analytical(a_val, T_inf_val, T_0_val, loud=False):
#    
#    # Define analytical symbols
#    a, t, T_inf, T_0 = sm.symbols("a t T_inf T_0")
#    T = sm.Function('T')
#    
#    # Define ode and intial value
#    ode = sm.Derivative(T(t),t) + a*(T(t)-T_inf)
#    
#    # Solve analytical
#    # Newer versions handle intial value problem more elegantly
#    if float(sm.__version__) >= 1.2:
#        ics = {T(0): T_0}
#        sol = sm.dsolve(ode, T(t), ics=ics )
#    else:    
#        sol = sm.dsolve(ode, T(t))
#        # Solve C1
#        C1 = sm.Symbol('C1')
#        C1_ic = sm.solve(sol.rhs.subs({t:0})-T_0,C1)[0]
#        sol   = sol.subs({C1:C1_ic})
#    
#    # Substitute constant numerical values
#    sub = sol.rhs.subs({a     : a_val,
#                        T_inf : T_inf_val,
#                        T_0   : T_0_val
#                        })
#    if loud:
#        print("ODE solution\t\t", sol)
#        print("After substitutions\t", sub)
#    
#    # Return function T(t)
#    # "numpy" adds 2 orders of magnitude np.ed
#    return sm.lambdify(t,sub,"numpy") 
#
#
##########################################################################
#if __name__ == "__main__":
#    
#    # paramteres
#    h       = 10
#    cp      = 910
#    rho     = 2700
#    d       = 2e-2
#    T_0     = 40
#    T_inf   = 20
#    
#    # Evaluation times
#    t = np.linspace(0, 60*60, 500)
#    
#    # Help vars 
#    V = 4/3*np.pi*(d/2)**3
#    A = 4*np.pi*(d/2)**2
#    m = rho*V
#    a = A*h/(m*cp)
#        
#    # Get analytical solution 
#    Tf = get_analytical(a, T_inf, T_0, loud=True)
#    
#    # Plot
#    plt.plot(t/60, Tf(t))
#    

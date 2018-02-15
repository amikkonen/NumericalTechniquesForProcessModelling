#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 14:26:16 2017

@author: ami
"""
import time
import scipy as sp
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI


# Globals
# usually a bad idea but good for truly global constants like gravity
g = 9.81


def Nu_func(k, nu, Pr, beta, dT, D):
    """
    Mills 1999, Basic Heat and Mass Transfer, page 294, Eq. 4.87
    
    Example 4.6 on page 300 for verification.
    """
    
    Ra = beta*dT*g*D**3/nu**2*Pr
    if Ra < 1e9 and Ra > 1e-6:
        Nu = 0.36 + 0.518*Ra**(1/4) / (1+(0.559/Pr)**(9/16))**(4/9)
    else:
        raise ValueError("Ra outside range. Ra=", Ra)
    return Nu

def h_func(Nu, k, D):
    return k/D*Nu

def q_func(h,  dT):
    return h*dT

    
############################################################################

def test_Nu_h_Q():
    """Verification of Nusselt number calculation.
    
    Mills 1999, Basic Heat and Mass Transfer, 
    
    Example 4.6 on page 300 for verification.
    
    """
    
    D      = 0.3        # m
    k      = 0.0331     # W/mK
    nu     = 25.5e-6    # m^2/s
    Pr     = 0.69
    beta   = 1/400
    dT     = 200        # K
    dx      = 1          # m
    A      = sp.pi*D*dx  # m^2
    
    Nu = Nu_func(k, nu, Pr, beta, dT, D)
    print()
    print("Nu", Nu)
    print("Should be 42.9")
    
    h = h_func(Nu, k, D)
    print()
    print("h", h)
    print("Should be 4.73")
    
    q = q_func(h, dT)
    Q = q*A
    print()
    print("Q", Q)
    print("Should be 892")
    
def main():
    # Pipe diameter (m)
    d      = 0.03      
    # Water intial temperature (C)
    T0     = 4
    # Outside temperature (C)
    T_inf  = -10

    
    # Fluid properies
#    k      = 0.0331      # W/mK
#    nu     = 25.5e-6     # m^2/s
#    Pr     = 0.69
#    rho    = 1000
#    cp     = 4.187e3    # J/kgK

    # Reference values
    p     = 1e5
    T_ref = 273 + (T0+T_inf) / 2
    
    k       = PropsSI('conductivity', "T", T_ref, "P", p, "air")
    mu      = PropsSI('V', "T", T_ref, "P", p, "air")
    rho_air = PropsSI('D', "T", T_ref, "P", p, "air")
    nu      = mu/rho_air
    Pr      = PropsSI('Prandtl', "T", T_ref, "P", p, "air")
    rho     = PropsSI('D', "T", T0+273, "P", p, "water")
    cp      = PropsSI('C', "T", T0+273, "P", p, "water")
    
    
    # Maximum time (s)
    seconds_in_hour = 60*60
    t_max  = 24*seconds_in_hour    
    
    # Time step (s)
    dt = 1
                      
    # Maximum number of steps 
    steps_max = round(t_max/dt)
    
    t_lst = []
    T_lst = []
    t_lst.append(0)
    T_lst.append(T0)
    
    t = 0
    # Loop through time
    for step in range(steps_max):
        t += dt
        
        # Help variables
        T_pipe   =  T_lst[-1]
        T_mean_K = 273 + (T_pipe + T_inf) / 2
        dT       = T_pipe - T_inf
        beta     = 1/T_mean_K
        
        # Heat transfer    
        Nu = Nu_func(k, nu, Pr, beta, dT, d)
        h  = h_func(Nu, k, d)
        q  = q_func(h,  dT)
        
        # Explicit Euler
        T_new = T_pipe - 4*q/(rho*d*cp)*dt
        
        # Check if new temperature below zero
        if T_new < 0:
            print("dt (s)", dt, "t (min)", t/60)
            break
        
        # Add new values to list
        T_lst.append(T_new)
        t_lst.append(t)
    T = sp.array(T_lst)

    # Plotting
    plt.plot(sp.array(t_lst)/60, T, label="dt="+str(dt))    
    plt.legend()    
    plt.axhline(0, color="k")
    plt.grid()
    plt.xlim(0,t/60)
    

############################################################################

if __name__ == "__main__":

    start = time.time()
    print("START")
#    test_Nu_h_Q()
    main()
    print("END", time.time() - start, "s")
    
    
    
    
    

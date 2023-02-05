#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of transient one-dimensional heat conduction.

KEB-45250
Numerical Techniques for Process Modelling,
Prosessien numeerinen mallinnus

Emphesis on simplisity. Low performance in larger systems.

Created on 2018 Jan 5 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg

##############################################################################

def solver(n, kt, L, T_B, q, T, rhoc, dt, t_max):
    
    # Cell size
    dx = L / n
    
    # Number of steps
    steps = int(np.ceil(t_max/dt))
    
    # Current time (s)
    t = 0
    
    # Coefficients
    # (Unnessasary to keep aW and aE appart, I just like it for generality)
    aW = kt/dx
    aE = kt/dx
    aP0 = rhoc*dx/dt
    
   # Source vector
    bConstant = np.zeros(n)
    # Coefficient matrix
    A = np.zeros((n,n))
    # Temperatures    
    Ts = np.zeros((steps+1,n))
    Ts[0,:] = T

    # Add aW to matrix A
    for k in range(1,n):
        A[k,k]   += aW
        A[k,k-1] += -aW
    
    # Add aE to matrix A
    for k in range(n-1):
        A[k,k]   += aE
        A[k,k+1] += -aE
    
    # Add time coefficient to diagonal
    for k in range(n):
        A[k,k] += aP0
    
    # Boundary A on the left
    A[-1,-1]      += kt/(dx/2)
    bConstant[-1] += kt/(dx/2)*T_B
    
    # Add heat generation
    bConstant += q*dx
    
    # Solution
    for step in range(1,steps+1):
        b  = bConstant + aP0*T
        T  = np.linalg.solve(A,b)    
        Ts[step] = T
        t += dt
    
    # Post
#    print("A\n", A)
#    print("b\n", b)
#    print("T\n",T)
    
    return Ts, dx
    
    
##############################################################################

if __name__ == "__main__":
    print("START")
    # Number of control volumes
    n = 5
    # Thermal conductivity (W/mK)
    kt = 10
    # Lenght (m)
    L = 0.02
    # Boundary temperatures (C)
    T_B = 0
    
    # Volumetric heat power (W/m3)
    q = 0#5e7

    # Initial temperatures (C)
    T_0 = 200

    # Product of density and heat capasity (J/m3K)
    rhoc = 10e6
    # Time step (s)
    dt = 2
    # Stop time (s)
    t_max = 1000
    # Plot at times (s)
    t_plot = [0,40,80,120]

    # Solve 
    Ts, dx = solver(n, kt,  L, T_B, q, T_0, rhoc, dt, t_max)    
    
    # Plot numerical and exact solution
    x = np.linspace(dx/2,L-dx/2,n)
    plt.figure(0)
    for t in t_plot:
        step = int(t/dt)
        plt.plot(x*1e3, Ts[step], "d-")
        print(str(t).ljust(4), Ts[step])
    plt.xlim(0,L*1e3)
    plt.xlabel("x (mm)")
    plt.ylabel("T (C)")
    
    # Plot first and last cell temperature
    plt.figure(1)
    plt.plot(np.arange(0,t_max+1,dt), Ts[:,0], label="first cell")
    plt.plot(np.arange(0,t_max+1,dt), Ts[:,-1], label="last cell")
    plt.legend()
    plt.xlim(0,t_max)
    plt.xlabel("t (s)")
    plt.ylabel("T (C)")    
    
    print("END")
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of transient one-dimensional heat conduction.

KEB-45250
Numerical Techniques for Process Modelling,
Prosessien numeerinen mallinnus

Exercice 4 problem 3

Emphesis on simplisity. Low performance in larger systems.

Created on 2018 Jan 5 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg

##############################################################################

def solver(n, kt,  L, T, T_A, T_B, crho, u, advection_scheme, q, dt, t_max):
    
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
    aP0 = crho*dx/dt
    
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
    
    #########################################################################
    # Advection
    ######################################################################### 
    
    # Central differencing
    if advection_scheme == "central_differencing":
        
        aWf = crho*u/2
        aEf = crho*u/2
        
        # Add aFw to matrix A
        for k in range(1,n):
            A[k,k]   -= aWf
            A[k,k-1] -= aWf
            
        # Add aFe to matrix A
        for k in range(n-1):
            A[k,k]   += aEf
            A[k,k+1] += aEf    
        
        # Boundary A on the left
        bConstant[0]   += crho*u*T_A
        
        # Boundary B on the right
        bConstant[-1]  += crho*u*T_B
    
    # Upwind
    elif advection_scheme == "upwind":

        aWf = crho*u
        aEf = crho*u
        
        # Positive velocity
        if u > 0:
            # Add aFw to matrix A
            for k in range(1,n):
                A[k,k-1] -= aWf
                
            # Add aFe to matrix A
            for k in range(n-1):
                A[k,k]   += aEf
            
            # Boundary A on the left
            bConstant[0]   += crho*u*T_A
            
            # Boundary B on the right
            A[-1,-1] += aEf        
        
        # Negative velocity
        elif u < 0:
            # Add aFw to matrix A
            for k in range(1,n):
                A[k,k] -= aWf
                
            # Add aFe to matrix A
            for k in range(n-1):
                A[k,k+1]   += aEf
            
            # Boundary A on the left
            A[0,0] -= aWf
            
            # Boundary B on the right
            bConstant[-1]   -= crho*u*T_B
            
        # No velocity    
        else:
            pass
    
    # Unkown advection scheme    
    else:
        raise ValueError("Unkown advection scheme")
    
    
    
    
    
    
    
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
    # Thermal conductivity (W/mK)
    kt = 0.1
    # Lenght (m)
    L = 1
    # Boundary temperatures (C)
    T_A = 10
    T_B = 5
    # Product of heat capacity and density
    crho = 1

#    advection_scheme = "central_differencing"
    advection_scheme = "upwind"
    
    # Velocity and Number of control volumes
    # Case a
    n = 50 ; u = 0.1
    # Case b
#    n = 5 ; u = 2.5
#     Case c
#    n = 20 ; u = 0
    
    # Case d
#    n = 5 ; u = -0.1
    # Case e
#    n = 5 ; u = -2.5
#     Case f
#    n = 20 ; u = -2.5
    
    # Volumetric heat power (W/m3)
    q = 0#5e7

    # Initial temperatures (C)
    T_0 = T_A

    # Time step (s)
    dt = 1
    # Stop time (s)
    t_max = 120
    # Plot at times (s)
    t_plot = [0,40,80,120]

    # Solve 
    Ts, dx = solver(n, kt,  L, T_0, T_A, T_B, crho, u, advection_scheme, 
                    q, dt, t_max)
    
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
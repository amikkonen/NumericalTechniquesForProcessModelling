#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of one-dimensional heat conduction.

KEB-45250
Numerical Techniques for Process Modelling,
Prosessien numeerinen mallinnus

Exercice 4 extra 2

Emphesis on simplisity. Low performance in larger systems.

Created on 2018 Jan 5 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy import linalg

##############################################################################

def solver(n, kt, Ac, L, T_A, T_B, q=0):
    
    # Cell size
    dx = L / n
    
    # Interpolate kt to faces
    ktf = (kt[1:] + kt[:-1]) / 2
    # First order extrapolation to boundaries
    ktfA = kt[0]
    ktfB = kt[-1]
    
    
    # Coefficients
    aW = ktf/dx*Ac
    aE = ktf/dx*Ac
    
    # Source vector
    b = sp.zeros(n)
    # Coefficient matrix
    A = sp.zeros((n,n))

    # Add aW to matrix A
    for k in range(1,n):
        A[k,k]   += aW[k-1]
        A[k,k-1] += -aW[k-1]

    # Add aE to matrix A
    for k in range(n-1):
        A[k,k]   += aE[k]
        A[k,k+1] += -aE[k]
    
    # Boundary A on the left
    A[0,0] += ktfA/(dx/2)*Ac
    b[0]   += ktfA/(dx/2)*Ac*T_A
    
    # Boundary B on the right
    A[-1,-1] += ktfB/(dx/2)*Ac
    b[-1]    += ktfB/(dx/2)*Ac*T_B
    
    # Add heat generation
    b += q*Ac*dx
    
    # Solution
    T = sp.linalg.solve(A,b)
    
    # Post
#    print("A\n", A)
#    print("b\n", b)
#    print("T\n",T)
    
    return T, dx
    
    
##############################################################################

if __name__ == "__main__":
    print("START")
    # Number of control volumes
    n = 20
    # Cross-sectional area (m2)
    Ac = 10e-3
    # Lenght (m)
    L = 0.5
    # Boundary temperatures (C)
    T_A = 100
    T_B = 500
    
    # Volumetric heat power (W/m3)
    q = 100 
    
    kt = sp.zeros(n)
    Told = sp.ones(n) * (T_A + T_B)/ 2
    # Iterative loop
    max_loops = 1000
    for k in range(max_loops):
        # Update kt        
        for kk in range(n):
            kt[kk] = PropsSI('conductivity','P',1e5,'T',Told[kk]+273.15,'Air')
        # Solve linear system    
        T, dx = solver(n, kt, Ac, L, T_A, T_B, q)
        # Check for convergence
        if sp.absolute((T - Told)).max() < 0.001:
            print("Converged after", k+1, "iterations")
            break
        Told = T
    
    # Check if loop terminated    
    if k == max_loops-1:
        print("WARNING! Maximum number of iterations reached!")
    
    # Plot numerical  solution
    plt.plot(sp.linspace(dx/2,L-dx/2,n), T, "k--o", label="numerical")
    plt.legend()
    plt.xlim(0,L)
    
    print("END")

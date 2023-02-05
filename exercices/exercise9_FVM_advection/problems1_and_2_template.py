#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of steady one-dimensional heat conduction and advection.

KEB-45250
Numerical Techniques for Process Modelling,
Prosessien numeerinen mallinnus

Emphesis on simplisity. Low performance in larger systems.

Created on 2018 Feb 9 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg

##############################################################################

def conduction(kt, dx, T_A, T_B, A, b):
    """Note that there is no need to return anything as Pyhton uses
       references for A and b. 
    """
    
    # Coefficients 
    aWc = kt/dx
    aEc = kt/dx
    
    # Add aW to matrix A
    for k in range(1,n):
        A[k,k]   += aWc
        A[k,k-1] += -aWc
    
    # Add aE to matrix A
    for k in range(n-1):
        A[k,k]   += aEc
        A[k,k+1] += -aEc
    
    # Boundary A on the left
    A[0,0] += kt/(dx/2)
    b[0]   += kt/(dx/2)*T_A
    
    # Boundary B on the right
    A[-1,-1] += kt/(dx/2)
    b[-1]    += kt/(dx/2)*T_B

def central_differencing(crho, u, T_A, T_B, A, b):
    """Note that there is no need to return anything as Pyhton uses
       references for A and b. 
    """

    aWf = crho*u/2
    aEf = crho*u/2

    # Add aWf to matrix A
    for k in range(1,n):
        A[k,k]   += -aWf
        A[k,k-1] += -aWf

    # Add aEf to matrix A
    for k in range(n-1):
        A[k,k]   += aEf
        A[k,k+1] += aEf

    # Boundary A on the left
    b[0]   += crho*u*T_A

    # Boundary B on the right
    b[-1]    += -crho*u*T_B






def upwind_scheme(crho, u, T_A, T_B, A, b):
    """Note that there is no need to return anything as Pyhton uses
       references for A and b. 
    """
    
    aWf = crho*u
    aEf = crho*u
    
    
    # Positive velocity
    if u > 0:

        # Add aWf to matrix A
        for k in range(1,n):
            A[k,k-1] += -aWf    

        # Add aEf to matrix A
        for k in range(n-1):
            A[k,k]   += aEf

        # Boundary A on the left
        b[0]   += crho*u*T_A

        # Boundary B on the right
        A[-1,-1]   += aEf


    # Negative velocity
    elif u < 0:
        # Add aWf to matrix A
        for k in range(1,n):
            A[k,k] += -aWf   
            
        # Add aEf to matrix A
        for k in range(n-1):
            A[k,k+1]   += aEf
            
        # Boundary A on the left
        A[0,0]   += -aWf
        
        # Boundary B on the right
        b[-1]   += -crho*u*T_B
        
    # No velocity    
    else:
        pass
    
def solver(n, kt,  L, T_A, T_B, crho, u, advection_scheme):
    
    # Cell size
    dx = L / n
  
    # Source vector
    b = np.zeros(n)
    # Coefficient matrix
    A = np.zeros((n,n))
    
    #########################################################################
    # Conduction
    #########################################################################    
    conduction(kt, dx, T_A, T_B, A, b)

    #########################################################################
    # Advection
    ######################################################################### 
    
    # Peclet number
    Pe = np.absolute(crho*u / (kt/dx))
    
    if advection_scheme == "central_differencing":
        central_differencing(crho, u, T_A, T_B, A, b)
    elif advection_scheme == "upwind":
        upwind_scheme(crho, u, T_A, T_B, A, b)
    # Unkown advection scheme    
    else:
        raise ValueError("Unkown advection scheme")
   
    #########################################################################
    # Solution
    ######################################################################### 
    
    # Solution
    T = np.linalg.solve(A,b)
    
    # Post
    print("Pe =", Pe)
    print("A\n", A)
    print("b\n", b)
    print("T\n",T)
    
    return T, dx
    
    
def analytical(n, kt, L, T_A, T_B, x, crho, u):
    return T_A + (T_B-T_A) * (np.exp(crho*u*x/kt)-1)/(np.exp(crho*u*L/kt)-1)
    
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

    advection_scheme = "central_differencing"
    # advection_scheme = "upwind"
    
    # Velocity and Number of control volumes
    # Case a
    # n = 5 ; u = 0.1
    # Case b
    n = 5 ; u = 2.5
    # Case c
    # n = 5 ; u = -0.1
    # Case d
    # n = 5 ; u = -2.5
    
    T, dx = solver(n, kt, L, T_A, T_B, crho, u, advection_scheme)
    
    x_exact = np.linspace(0,L)
    T_exact = analytical(n, kt, L, T_A, T_B, x_exact, crho, u)
    
    # Plot numerical and exact solution
    x = np.linspace(dx/2,L-dx/2,n)    
    plt.plot(x_exact, T_exact, "r", label="exact")
    plt.plot(x, T, "k--o", label="numerical")
    plt.legend()
    plt.xlim(0,L)
    
    print("END")
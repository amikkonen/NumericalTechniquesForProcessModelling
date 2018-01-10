#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of steady one-dimensional heat conduction.

See Fersteeg and Malalasekera 2007, example 4.1 om page 118 for details.

No sources, constant propersties, contant temperature boundaries.

Emphasis on simplicity. Low performance in larger systems.

Created on Thu Jan  4 15:48:25 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from scipy import linalg
from matplotlib import pyplot as plt

##############################################################################

def main(n, kt, Ac, L, T_A, T_B):
    
    # Cell size
    dx = L / n
    # Coefficients
    aW = kt/dx*Ac
    aE = kt/dx*Ac
    
    # Source vector
    b = sp.zeros(n)
    # Coefficient matrix
    A = sp.zeros((n,n))

    # Add aW to matrix A
    for k in range(1,n):
        A[k,k]   += aW
        A[k,k-1] += -aW
#    print("A with aW\n", A)
    
    # Add aE to matrix A
    for k in range(n-1):
        A[k,k]   += aE
        A[k,k+1] += -aE
#    print("A with aE\n", A)
    
    # Boundary A on the left
    A[0,0] += kt/(dx/2)*Ac
    b[0]   += kt/(dx/2)*Ac*T_A
    
    # Boundary B on the right
    A[-1,-1] += kt/(dx/2)*Ac
    b[-1]    += kt/(dx/2)*Ac*T_B
    
    # Solution
    T = sp.linalg.solve(A,b)
    
    # Post
    print("A\n", A)
    print("b\n", b)
    print("T\n",T)
    
    # Plot numerical and exact solution
    plt.plot([0,L], [T_A, T_B], "r", label="exact")
    plt.plot(sp.linspace(dx/2,L-dx/2,n), T, "k--o", label="numerical")
    plt.legend()
    plt.xlim(0,L)
    
##############################################################################

if __name__ == "__main__":
    print("START")
    # Number of control volumes
    n = 5
    # Thermal conductivity (W/mK)
    kt = 1000
    # Cross-sectional area (m2)
    Ac = 10e-3
    # Lenght (m)
    L = 0.5
    # Boundary temperatures (C)
    T_A = 100
    T_B = 500
    
    main(n, kt, Ac, L, T_A, T_B)
    print("END")

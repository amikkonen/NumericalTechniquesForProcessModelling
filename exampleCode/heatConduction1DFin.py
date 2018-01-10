#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of steady one-dimensional heat conduction in a fin.

See Fersteeg and Malalasekera 2007, example 4.3 om page 125 for details.

Linear source (heat convection), constant propersties, 
contant temperature boundary, zero gradient (insulated) boundary.

Emphasis on simplicity. Low performance in larger systems.

Created on Thu Jan  4 16:54:30 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from scipy import linalg
from matplotlib import pyplot as plt

##############################################################################

def fin(L, T_B, T_inf, n2, n):
    # Cell size
    dx = L / n
    
    # Coefficients
    aW = 1/dx
    aE = 1/dx
    
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
    
    # Add linear source to matrix 
    for k in range(n):
        A[k,k] += n2*dx
    
    # Add linear source to vector b
    for k in range(n):
        b[k] += n2*dx*T_inf
        
    # Boundary A on the left
    A[0,0] += 2/dx
    b[0]   += 2/dx*T_B

    # Solution
    T = sp.linalg.solve(A,b)
    
#    print("A\n", A)
#    print("b\n", b)
#    print("T\n",T)
    
    return T, dx

##############################################################################

if __name__ == "__main__":
    print("START")
    # Number of control volumes
    ns = [5,50]
    # Lenght (m)
    L = 1
    # Boundary temperatures (C)
    T_B   = 100
    # Enviroment temperature
    T_inf = 20
    # Fin paramter (1/m2)
    n2 = 25
    
    T  = []
    dx = []
    for k, n in enumerate(ns):
        Ti, dxi = fin(L, T_B, T_inf, n2, n)
        T.append(Ti)
        dx.append(dxi)

    # Exact solution
    x = sp.linspace(0,L)
    T_exact = T_inf + (T_B-T_inf)*sp.cosh(n2**0.5*(L-x))/sp.cosh(n2**0.5*L) 
    
    #  Plot numerical and exact solution
    plt.plot(x, T_exact, "k", label="exact")
    plt.plot(sp.linspace(dx[1]/2,L-dx[1]/2,ns[1]), T[1], "d", label="n="+str(ns[1]))    
    plt.plot(sp.linspace(dx[0]/2,L-dx[0]/2,ns[0]), T[0], "o", label="n="+str(ns[0]))
    plt.legend()
    plt.xlim(0,L)
    print("END")

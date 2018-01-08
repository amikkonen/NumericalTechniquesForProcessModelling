#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A course example of unsteady one-dimensional heat conduction.

See Fersteeg and Malalasekera 2007, example 8.2 om page 253 for details.

Constant propersties, contant temperature boundary, zero gradient (insulated) 
boundary, transient.

Emphasis on simplicity. Low performance in larger systems.

Created on Thu Jan  4 17:45:12 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import scipy as sp
from matplotlib import pyplot as plt

##############################################################################

def solver(L, T_B, T, ks, rhoc, n, dt, t_max):
    # Cell size
    dx = L / n
    
    # Current time (s)
    t = 0
    
    # Number of steps
    steps = int(sp.ceil(t_max/dt))
    
    # Coefficients
    aW  = ks/dx
    aE  = ks/dx
    aP0 = rhoc*dx/dt
    
    # Source vector
    bConstant = sp.zeros(n)
    # Coefficient matrix
    A = sp.zeros((n,n))
    # Temperatures    
    Ts = sp.zeros((steps+1,n))
    Ts[0,:] = T
    
    # Add aW to matrix A
    for k in range(1,n):
        A[k,k]   += aW
        A[k,k-1] += -aW
    
    # Add aE to matrix A
    for k in range(n-1):
        A[k,k]   += aE
        A[k,k+1] += -aE
    
    # Constant T boundary on the right
    A[-1,-1]      += 2*ks/dx
    bConstant[-1] += 2*ks*T_B

    # Add time coefficient to diagonal
    for k in range(n):
        A[k,k] += aP0

    # Solution
    for step in range(1,steps+1):
        b  = bConstant + aP0*T
        T  = sp.linalg.solve(A,b)    
        Ts[step] = T
        t += dt
    
#    print("A\n", A)
#    print("bConstant\n", bConstant)
#    print("b\n", b)
#    print("T\n",T)
    
    return Ts, dx

def exact_solution(L, T_B, T_0, ks, rhoc, n, t, x):
    """Özisik (1985)
    
    NOTE: This correlation doesn't have T_B as a parameter. It is assumed here
          that the correlation only works if T_B == 0. Not checked.
    
    SOMETHING WRONG!
    
    """
    # Check if T_B == 0
    assert T_B == 0 
    
    # Sum terms
    n = sp.arange(1,n)
    
    lambdan = (2*n-1)*sp.pi/(2*L)
    alpha   = ks/rhoc

    T = sp.zeros(len(x))
    for k in range(len(x)):
        T[k] = T_0 + 4/sp.pi * ((-1)**(n+1)/(2*n-1)*\
                         sp.exp(-alpha*lambdan**2*t)*sp.cos(lambdan*x[k])).sum()
    return T 
    

##############################################################################

if __name__ == "__main__":
    print("START")

    # Number of control volumes
    n = 5
    # Lenght (m)
    L = 0.02
    # Boundary temperatures (C)
    T_B   = 0
    # Initial temperatures (C)
    T_0 = 200
    # Thermal conductivity (W/mK)
    ks = 10 
    # Product of density and heat capasity (J/m3K)
    rhoc = 10e6
    # Time step (s)
    dt = 2
    # Stop time (s)
    t_max = 120
    # Plot at times (s)
    t_plot = [0,40,80,120]

    Ts, dx = solver(L, T_B, T_0, ks, rhoc, n, dt, t_max)
    
#    x_exact = sp.linspace(dx/2,L-dx/2,50)
    x       = sp.linspace(dx/2,L-dx/2,n)
    for t in t_plot:
        step = int(t/dt)
#        T_exact = exact_solution(L, T_B, T_0, ks, rhoc, 10000, t, x_exact)
#        plt.plot(x_exact, T_exact, "k-")
        plt.plot(x, Ts[step], "d")
        print(str(t).ljust(4), Ts[step])

    print("END")

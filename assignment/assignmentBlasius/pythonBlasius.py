#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:05:27 2018

@author: ami
"""
import scipy as sp
from scipy.integrate import solve_bvp
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

def dimmless_blasius(eta_max=10, n=1000):
    
    def odefun(eta, F):
        f1 = F[0]
        f2 = F[1]
        f3 = F[2]
        return sp.vstack([f2,
                         f3,
                         -0.5 * f1 * f3])

    def bcfun(fa, fb):
        return sp.array([fa[0],      # f1(0) =  0
                        fa[1],       # f2(0) = 0
                        1.0 - fb[1]  # f2(inf) = 1
                        ]) 

    eta = sp.linspace(0, eta_max, n)
    Finit = sp.vstack([eta, sp.exp(-eta), sp.exp(-eta)])

    sol = solve_bvp(odefun, bcfun, eta, Finit)
    
#    f   = sol.y[0]
    fp  = sol.y[1]
#    fpp = sol.y[2]
    
    # u = fp*U
    return eta, fp
    

def blasius_plate(U, x, y, nu, eta_max=7, n=1000):
    eta_1d, fp_1d = dimmless_blasius(eta_max, n)
    fp_interp =  interp1d(eta_1d, fp_1d,
                          bounds_error=False,
                          fill_value=(fp_1d[0],fp_1d[-1]))
    X, Y = sp.meshgrid(x,y)
    eta = Y * sp.sqrt(U/(nu*X))
    fp = fp_interp(eta)
    u = U*fp
    return u    

if __name__ == "__main__":
    print("START")
    U = 0.2
    x = sp.linspace(1e-20,0.1,5)
    y = sp.linspace(0,0.02,1000)
    nu = 15.11e-6 # air at 20C
    
    
    eta_max=7
    n = 1000
    
    Re = U*x.max() / nu
    print("Re", Re)
    
    u = blasius_plate(U, x, y, nu, eta_max, n)
    
    for k in range(len(x)):
        plt.plot(u[:,k],y, '-')
    
    print("END")


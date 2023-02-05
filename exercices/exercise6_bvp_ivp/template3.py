#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:05:27 2018

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
import numpy as np
from scipy.integrate import solve_bvp
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

def dimmless_blasius(eta_max=10, n=1000):
    
    def ode(eta, y):
        return np.array([,
                         ,
                         )

    def bc(ya, yb):
        return np.array([, 
                        ,  
                             
                        ]) 

    eta   = np.linspace(0, eta_max, n)
    yinit = np.vstack([eta, np.exp(-eta), np.exp(-eta)])

    sol = solve_bvp(ode, bc, eta, yinit)
    
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
    X, Y = np.meshgrid(x,y)
    eta = Y * np.sqrt(U/(nu*X))
    fp = fp_interp(eta)
    u = U*fp
    return u    

if __name__ == "__main__":
    print("START")
    U = 0.2
    x = np.linspace(1e-20,0.1,5)
    y = np.linspace(0,0.02,1000)
    nu = 15.11e-6 # air at 20C
    
    
    Re = U*x.max() / nu
    print("Re_x max", Re)
    
    u = blasius_plate(U, x, y, nu)
    
    for k in range(len(x)):
        plt.plot(u[:,k],y*1e3, '-', label="x="+str(int(x[k]*1e3))+"mm")
    
    print(u[3,3])
    
    plt.xlabel("U (m/s)")
    plt.ylabel("y (mm)")
    plt.legend(frameon=False)
    print("END")


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 13:59:39 2018

@author: ami
"""
import scipy as sp
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

def main():
    g = 9.81
    
    def free_fall_analytical(t, y0=0, v0=0):
        return sp.array([y0+v0*t-0.5*g*t**2,
                         v0+g*t,
                        ])
    
    def free_fall(t, y):
        return sp.array([y[1], -9.8])

    # Last time    
    t_end = 5
    
    # Initial position and velocity
    xv_inti = [0,10]
    
    # Times to solve the analytical system
    t_ana = sp.linspace(0, t_end, 100)
    
    # Times to visualize the numerical solution. Nothing to do with timestep.
    t_eval = [0,2,4,t_end]
    
    # Numerical solution
    sol = solve_ivp(free_fall, [0, t_end], xv_inti,
                    t_eval=t_eval
                    )
    
    # Post
    print("\nt\n",sol.t)
    print("\ny\n", sol.y[0])
    print("\nV\n", sol.y[1])
    
    # Plot
    plt.plot(t_ana, free_fall_analytical(t_ana,xv_inti[0], xv_inti[1])[0])
    plt.plot(sol.t, sol.y[0], 'd')

if __name__ == "__main__":
    print("START")
    main()
    print("END")

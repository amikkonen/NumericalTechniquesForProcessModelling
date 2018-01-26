#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:09:38 2018

@author: ami
"""

import scipy as sp
from scipy import linalg
from matplotlib import pyplot as plt
import sympy

#from sympy import solve, symbols, pi, Eq
#from sympy.physics.units import Quantity, length, mass
#from sympy.physics.units import day, gravitational_constant as G
#from sympy.physics.units import meter, kilogram
#
#
#T = symbols("T")
#a = Quantity("a", length, 108208000e3*meter)
#M = Quantity("M", mass, 1.9891e30*kilogram)
#eq = Eq(T**2 / a**3, 4*pi**2 / G / M)
#q = solve(eq, T)
#
#print(q)

#import sympy as sm
#x,t,z,nu = sm.symbols('x t z nu')
##sm.init_printing(use_unicode=True)
##print(sm.diff(sm.sin(x)*sm.exp(x),x))
#
#
#y = sm.Function('y')
#print(sm.dsolve(sm.Eq(y(t).diff(t, t) - y(t), sm.exp(t)), y(t)))
#print(sm.latex(sm.Integral(sm.cos(x)**2, (x, 0, sm.pi))))

#6.28318530717959*venus_a**(3/2)/(sqrt(gravitational_constant)*sqrt(solar_mass))

#x = sp.linspace(-10,10)
#
#a = 15
#b = -40
#y = a*x + b
#plt.plot(x,y,'--', label="y=15x-40")
#plt.axhline(0, color='k')
#
#a = 2
#b = -3
#c = -30
#
#y2 = a*x**2 + b*x + c 
#plt.plot(x,y2, label="$y=2x^2-3x-30$")
#
#
#
#y3 = sp.sin(x/sp.pi)*100+sp.sinh(x/sp.pi)*15
#plt.plot(x,y3, 'd',label="$100sin(x/pi)+15sinh(x/pi)$")
#
#
#
#
#plt.legend()
#plt.grid()






#def f(Re, d, eps=0):
 
Re  = 10000
d   = 1e-2
eps = 0
    
f = 1e-2    
n = 10
for k in range(n):
    print(f)
    f = (-2*sp.log10(eps/3.7/d
         +2.51/(Re*f**0.5))
        )**-2     
    

#A = sp.array([[3,2,-1],
#              [2,-2, 4],
#              [-1,0.5,-1],
#              ])
#b = sp.array([1,-2,0])    
#x = linalg.solve(A,b)
#print(A@x-b)


#
#
#b = sp.array([1,-2,0])    
#x = sp.array([2,1,1])
#
#n = 100
#for k in range(n):
#    A = sp.array([[3*x[2],2,-1],
#                  [2,-2, 4],
#                  [-1,0.5,-1],
#                  ])
#    x = linalg.solve(A,b)
#
## Print solution
#print(x)
## Test solution, @ is matrix multiplication
#print(A@x-b)    
#
#
##b = sp.array([1,-2,0])    
##x = sp.array([2,1,1])
##
##relax = 0.5
##n = 100
##for k in range(n):
##    A = sp.array([[3*x[1],2,-1],
##                  [2,-2, 4],
##                  [-1,0.5,-1],
##                  ])
##    x = relax*linalg.solve(A,b) + (1-relax)*x
##
##print(x)
##print(A@x-b)    
#    
#    
#
#
#
#def Ra_func(beta, dT, L, nu, alpha):
#    return beta*dT*9.81*L**3/(nu*alpha)


    
    
    
    
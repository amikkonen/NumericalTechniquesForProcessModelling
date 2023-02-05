#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:44:42 2021

@author: ami
"""

import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

def func(x):
    return x*np.sin(x) + 0.2*x**2 - 3

def negfunc(x):
    return -1*func(x)

def hfunc(x):
    return func(x+6)**0.5 - x 

def absfunc(x):
    return np.absolute(func(x))




def penalty_func(x):
    f = func(x)
    h = hfunc(x)
    
    if h < 0:
        f += 100
        
    # mask = h < 0
    # # print(mask)
    # f[mask] += 10
        
    return f
    


# 
x = np.linspace(0, 10,100)

plt.plot(x, func(x))
# plt.plot(x, negfunc(x), "--")

# plt.plot(x, penalty_func(x))
plt.plot(x, absfunc(x))


plt.plot(x, hfunc(x), "r-")
plt.axhline(0, color="k")
plt.xlabel("x")
plt.ylabel("func(x)")

#########################################################################
# (a) 
#########################################################################

# res = optimize.root(func, 5)
# print(res)
# assert(res.success)
# # f = res.x[0]

# plt.plot(res.x, res.fun, "d", label="root")



#########################################################################
# # (b) 
# #########################################################################

# res = optimize.root_scalar(func, bracket=[0, 10], method='brentq')
# print(res)
# assert(res.converged)

# plt.plot(res.root, 0, "x", label="root_scalar")



# #########################################################################
# # (c) 
# #########################################################################

# res = optimize.minimize_scalar(func, bounds=(0, 10), method='bounded')
# print(res)
# assert(res.success)

# plt.plot(res.x, res.fun, ".", label="minimize_scalar")


# #########################################################################
# # (d) 
# #########################################################################

# res = optimize.minimize_scalar(negfunc, bounds=(0, 10), method='bounded')
# print(res)
# assert(res.success)

# plt.plot(res.x, -res.fun, "o", label="max")

# # #########################################################################
# # # (f) 
# # #########################################################################

# res = optimize.minimize_scalar(penalty_func, bounds=(0, 10), method='bounded')
# print(res)
# assert(res.success)

# plt.plot(res.x, res.fun, "v", label="penalty")

# # print(func(0))

# #########################################################################
# # (e) 
# #########################################################################

res = optimize.minimize_scalar(absfunc, bounds=(3, 10), method='bounded')
print(res)
assert(res.success)

plt.plot(res.x, res.fun, "^", label="root with min")

# print(func(0))

###########################################################################



plt.legend(frameon=False)













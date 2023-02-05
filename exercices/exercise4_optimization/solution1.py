# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

def func(x):
    return x*np.sin(x) + 0.2*x**2 - 3

def negative_func(x):
    return -1*func(x)

def absolute_func(x):
    return np.absolute(func(x))

def hfunc(x):
    return func(x+6)**(1/2) - x

def penalty_func(x):
    f = func(x)
    h = hfunc(x)
    
    # Penalty
    if h < 0:
        f += 100
        
    return f
    
############################################################################
# (a) Root with root
############################################################################

res = optimize.root(func, 5)
#print(res)
assert(res.success)
plt.plot(res.x, res.fun, "o", label="root with root",markersize=12)

############################################################################
# (b) Root with root_scalar using bracket
############################################################################

res = optimize.root_scalar(func, bracket=[0,10])
# print(res)
assert(res.converged)
plt.plot(res.root, 0, "x", label="root with root_scalar", markersize=12)

############################################################################
# (c) Minimize with minimize_scalar
############################################################################

res = optimize.minimize_scalar(func, bounds=[0,10], method="bounded")
# print(res)
assert(res.success)
print("min x", res.x)
print("min f", res.fun)
plt.plot(res.x, res.fun, "d", label="minimize")

############################################################################
# (d) Maximize with minimize_scalar
############################################################################

res = optimize.minimize_scalar(negative_func, bounds=[0,10], method="bounded")
# print(res)
assert(res.success)
plt.plot(res.x, -res.fun, "d", label="maximize")

############################################################################
# (e) Root with minimize_scalar
############################################################################

res = optimize.minimize_scalar(absolute_func, bounds=[3,10], method="bounded")
# print(res)
assert(res.success)
plt.plot(res.x, res.fun, "*", label="root with minimize", markersize=12)

############################################################################
# (f) Minimize with minimize_scalar and constraint
############################################################################

res = optimize.minimize_scalar(penalty_func, bounds=[0,10], method="bounded")
# print(res)
assert(res.success)
print("(f) Minimize with minimize_scalar and constraint")
print("min x", res.x)
print("min f", res.fun)
plt.plot(res.x, res.fun, "d", label="minimize with constraint")

############################################################################
# (extra) Pyswarm
# Will not run if pyswarm is not installed
############################################################################
try:
    from pyswarm import pso
    
    print("pso")
    # No constraint
    xopt, fopt = pso(func, 0, 10)
    
    # Constraint
    # xopt, fopt = pso(func, 0, 10, f_ieqcons=hfunc)
    
    print("min x", xopt)
    print("min f", fopt)
    
    plt.plot(xopt, fopt, "v", label="minimize pso")
    
except:
    print("No pyswarm available. Should be easy to install with pip.\nSee https://pythonhosted.org/pyswarm/")


############################################################################
# Plotting the full range
############################################################################
x = np.linspace(0,10)
plt.plot(x,func(x),"b",label="func")
hval = hfunc(x)
mask = hval>=0
plt.plot(x,hval,"r",label="h constraint")
plt.fill_between(x[mask],hval[mask],color="r",alpha=0.5)

# Prettyfy plot
plt.xlim(0,10)
plt.axhline(0, color="k")
plt.legend(frameon=False,prop={'size': 11})

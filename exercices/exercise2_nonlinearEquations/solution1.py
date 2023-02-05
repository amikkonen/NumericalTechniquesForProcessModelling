# -*- coding: utf-8 -*-
# import scipy as sp
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

def haaland(Re, eps_over_d):
    """ Friction factor, Haaland equation
    """
    f  = (-1.8*np.log10((eps_over_d/3.7)**1.11+6.9/Re))**-2
    return f

def naive_white_colebrook(Re, eps_over_d):
    """Simple iteration.
    
    """
    f_old = np.inf
    f = 1e-2
    while np.absolute(f-f_old) > 1e-8:
#        print(f)
        f_old = f
        f = (-2*np.log10( eps_over_d/3.7 + 2.51/(Re*f**0.5) ))**-2
    return f

def white_colebrook(Re, eps_over_d):
    """Friction factor, White-Colebrook equation
        Valid for Re > 4000 and Re < 5e6
    
    """

    def to_zero(f):
#        print(f)
        return f**(-0.5) + 2*np.log10( eps_over_d/3.7 + 2.51/(Re*f**0.5) )
    
#    # root, needs intial guess, gradient based, potentially unstable
#    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
    f_initial = 1e-2
    res = optimize.root(to_zero, f_initial)
#    print(res)
    assert(res.success)
    f = res.x[0]

    # root_scalar is a preferred solver for scalar valued functions.
    # "brentq" is a very reliable method.
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html
    # res = optimize.root_scalar(to_zero, bracket=[1e-12, 0.1], method="brentq")
#    print(res)
    # assert(res.converged)
    # f = res.root
    
    return f

##############################################################################
# (a) Naive iterations single value
##############################################################################
#Re = 1e4
#eps_over_d = 0
#print("Haaland f", haaland(Re, eps_over_d))
#print("White-Colebrook f", naive_white_colebrook(Re, eps_over_d))

##############################################################################
# (b) Single value
##############################################################################
#Re = 1e4
#eps_over_d = 0
#print("Haaland f", haaland(Re, eps_over_d))
#print("White-Colebrook f", white_colebrook(Re, eps_over_d))

##############################################################################
# (c) Re array
##############################################################################
# Reynolds numbers
# https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.linspace.html    
Re = np.linspace(4000, 5e6, 1000)
# Wall roughness
eps_over_d = 0

# Friction factors, initilized with as many zeros as there are Re values
f = np.zeros_like(Re)
# Loop over Re values
for k in range(len(Re)):
    f[k] = white_colebrook(Re[k], eps_over_d)
plt.loglog(Re, f, label="White-Colebrook")
plt.loglog(Re, haaland(Re, eps_over_d), label="Haaland")
plt.legend()


# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 11:42:50 2019

@author: niemel34
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

def K(T):
    """Equilibrium constant for water-gas-shift reaction CO + H2O = CO2 + H2O
       Taken from www.energy.kth.se, source Twigg, 1989
    """    
    z = 1000/T - 1
    return np.exp(0.31688 + 4.1778*z + 0.63508*z**2 - 0.29353*z**3)
    
    
def equilibrium_composition(a, T):
    """Function that solves composition from equilibrium and balance equations
    
    """
    
    # Given
    n_N2 = 2*a*3.76
    
    def to_zero(n):
        f = np.zeros(4)
        f[0] = n[0] + n[1] - 1
        f[1] = 2*n[0] + n[1] + n[2] - 4*a
        f[2] = n[2] + n[3] - 2
        f[3] = K(T)*n[1]*n[2] - n[0]*n[3]
        return f
    
    # Solve n
    # n[0] = n_CO2, n[1] = n_CO, n[2] = n_H2O, n[3] = n_H2
    sol = optimize.root(to_zero, 0.5*np.ones(4))
    
    assert sol.success
    n = sol.x
    
    # Add given n_N2 to n[4]
    n = np.append(n, n_N2)
    return n

##############################################################################
# (a) a = λ = 0.8, T = 1600K
##############################################################################
# # Air-fuel ratio a = lambda, value 1 mean stoichiometric composition
# a = 0.8
# # Temperature where combustion products are, unit Kelvin
# T = 1600 
# # Moles of each gas
# # n[0] = n_CO2, n[1] = n_CO, n[2] = n_H2O, n[3] = n_H2, n[4] = n_N2
# n = equilibrium_composition(a, T)
# # Mole percentage for each gas (and volume percentage if assumed ideal gas) 
# vol_percentage = 100 * ( n/sum(n) )

# # Out
# print('----------------------------------------------------------------')
# print('Equilibrium composition is:')
# print('\t', 'CO2:', round(vol_percentage[0],1), 'vol-%')
# print('\t', 'CO: ', round(vol_percentage[1],1), 'vol-%')
# print('\t', 'H2O:', round(vol_percentage[2],1), 'vol-%')
# print('\t', 'H2: ', round(vol_percentage[3],1), 'vol-%')
# print('\t', 'N2: ', round(vol_percentage[4],1), 'vol-%')
# print('----------------------------------------------------------------')


##############################################################################
# (b) 0.3 ≤ a ≤ 1, T = 1600K
##############################################################################

# Air-fuel ratio a = lambda, value 1 mean stoichiometric composition
a = np.linspace(0.3, 1, 50)
# Temperature where combustion products are, unit Kelvin
T = 1600 
# Moles of each gas
# n[0] = n_CO2, n[1] = n_CO, n[2] = n_H2O, n[3] = n_H2, n[4] = n_N2
# 2D array, 5 rows, len(a) cols
n = np.zeros([5, len(a)])
for k in range(len(a)):
    # : means all values, i.e. all rows here
    # equilibrium_composition returns an array with len 5
    # k loops through all air-fuel ratio cases (index)
    n[:, k] = equilibrium_composition(a[k], T)

# Mole percentage for each gas (and volume percentage if assumed ideal gas) 
# Without additional info, rum operates in "col" direction
vol_percentage = 100 * ( n/sum(n) )

# Plot results
names = ['CO2', 'CO', 'H2O', 'H2', 'N2']
for k in range(5):
    plt.plot(a, vol_percentage[k], label=names[k])

# Prettyfy    
plt.xlabel('Equivalence ratio λ (-)')
plt.ylabel('Volume percentage (%)')
plt.legend()

      

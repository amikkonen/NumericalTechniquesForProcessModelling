# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 11:50:58 2019

@author: Niko Niemelä, modified by Antti Mikkonen
"""
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI

C2K = 273.15

def h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid):
    """Zukauskas correlation for inline tube banks
    
    """
    
    # Constants for Re = 1e3 - 2e5
    C = 0.27 
    m = 0.63
    
    # Material properties in fluid mean temperature T_m
    T_m = (T_in + T_out) / 2
    mu = PropsSI('V', 'T', T_m, 'P', p, fluid)
    rho = PropsSI('D', 'T', T_m, 'P', p, fluid)
    nu = mu/rho
    Pr  = PropsSI('Prandtl', 'T', T_m, 'P', p, fluid)
    Pr_s  = PropsSI('Prandtl', 'T', T_s, 'P', p, fluid)
    k = PropsSI('conductivity', 'T', T_m, 'P', p, fluid)
    
    # Calculate Nusselt number and return h
    U_max = S_T*U_inf / (S_T - D)
    Re = U_max*D/nu
    if Re < 1e3 or Re > 2e5:
        print('Reynolds number not valid for correlation')
    Nu_D = C * Re**m * Pr**0.36 * (Pr/Pr_s)**0.25
    return Nu_D*k / D


def dT_log(T_in, T_out, T_s):
    """Logartimic temperature difference
    
    """
    dT_log = ( (T_s-T_in)-(T_s-T_out) ) / np.log( (T_s-T_in)/(T_s-T_out) )
    return dT_log


def iterate_solution(T_out, max_iter, N, S_T, U_inf, 
                            D, T_s, T_in, p, fluid, T_history=None):
    """Solve using under relaxation
    
    Retun heat power per unit legth of heat exhanger and outlet temperature 
    
    """
    
    rho_in = PropsSI('D', 'T', T_in, 'P', p, fluid)
    m = rho_in*U_inf*N_T*S_T
    # Under Relaxation Factor
    URF = 0.1
    
    for i in range(max_iter):
        # Heat capasity
        c_p = PropsSI('C', 'T', (T_in+T_out)/2, 'P', p, fluid)
        # Save the old T_out from last iteration round
        T_out_old = T_out 
        # Heat trasfer coefficient
        h = h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid) 
        # Heat power per unit legth
        Q = N*h*np.pi*D*dT_log(T_in, T_out, T_s)
        # Calculate new T_out
        T_out = T_in + Q/(m*c_p)
        # Calculate error and stop the for-loop if within tolerance
        change = abs( T_out - T_out_old )
        if change < 1e-7: # Limit for difference between two iterations (K)
            print('Iterations required:', i)
            break
        # Use underrelaxation for next iteration round
        T_out = URF*T_out + (1-URF)*T_out_old
    
        # Extra paramter for iteration history
        if T_history is not None:
            print("T_out%i"%(i+1), T_out)
            T_history.append(T_out)
    
        # Sanity check, temperature must be real 
        assert np.isreal(T_out), "T_out=" + str(T_out)

    # Make sure that you notice if the iteration did not converge
    # Will raise an AssertionError if i < max_iter - 1 is false
    assert i < max_iter - 1, "Did not converge"
    
    return T_out


def calculate_solution_with_scipy(N, S_T, U_inf, D, T_in, T_s, p, fluid):
    """Solve using scipy
    
    Return outlet temperature 
    
    """
    rho_in = PropsSI('D', 'T', T_in, 'P', p, fluid)
    m = rho_in*U_inf*N_T*S_T

    def energy_balance(T_out):
        T_m   = (T_in+T_out)/2
        c_p   = PropsSI('C', 'T', T_m, 'P', p, fluid)
        h     = h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid) 
        Q     = N*h*np.pi*D*dT_log(T_in, T_out, T_s)
        inbalance = Q - m*c_p*(T_out-T_in)
        return inbalance

    # Different calls
    T_out = optimize.brentq(energy_balance, T_s+0.1, T_in-0.1)
    return T_out


##############################################################################
# MAIN
##############################################################################


# Initial values for calculations
D       = 0.0635
S_T     = 2.5*D
N_T     = 24
N_L     = 100
T_in    = 900
T_s     = 600
U_inf   = 5
p       = 1e5
fluid   = 'Air'

N = N_T*N_L

# (a) With under-relaxation   
# Extra parameter for iteration history
T_history = []

max_iter = 500 # Maximum number of iterations
T_out0   = (T_s+T_in)/2 # Initial guess for outlet temperature
T_out = iterate_solution(T_out0, max_iter, N, S_T, U_inf, 
                            D, T_s, T_in, p, fluid, T_history)

plt.plot(T_history)

# (b) With scipy or similar
# T_out = calculate_solution_with_scipy(N, S_T, U_inf, D, T_in, T_s, p, fluid)

# Print results
# The following syntax means that %.1f is replaced T_out using 1 decimal place
print("T_out %.1f K" % T_out)



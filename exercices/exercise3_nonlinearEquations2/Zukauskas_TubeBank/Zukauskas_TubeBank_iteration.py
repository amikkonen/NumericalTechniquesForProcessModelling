# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 11:50:58 2019

@author: Niko Niemelä
"""
import scipy as sp
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI

#----  Zukauskas correlation for inline tube banks   ----

def h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid):
    
    # Constants for Re = 1e3 - 2e5
    C = 0.27 
    m = 0.63
    
    # Material properties in fluid mean temperature T_m
    T_m = (T_in+T_out)/2
    mu = PropsSI('V', 'T', T_m, 'P', p, fluid)
    rho = PropsSI('D', 'T', T_m, 'P', p, fluid)
    nu = mu/rho
    Pr  = PropsSI('Prandtl', 'T', T_m, 'P', p, fluid)
    Pr_s  = PropsSI('Prandtl', 'T', T_s, 'P', p, fluid)
    k = PropsSI('conductivity', 'T', T_m, 'P', p, fluid)
    
    # Calculate Nusselt number and return h
    U_max = S_T*U_inf / (S_T-D)
    Re = U_max*D/nu
    if Re < 1e3 or Re > 2e5:
        print('Reynolds number not valid for correlation')
    Nu_D = C * Re**m * Pr**0.36 * (Pr/Pr_s)**0.25
    return Nu_D*k / D

#---------------------  dT_log   --------------------------

def dT_log(T_in, T_out, T_s):
    
    dT_log = ( (T_s-T_in)-(T_s-T_out) ) / sp.log( (T_s-T_in)/(T_s-T_out) )
    return dT_log

#---------------------  Iterate solution   --------------------------
def iterate_solution(T_out0, N_max, N, S_T, U_inf, D, T_s, T_in, p, fluid):
    
    rho_in = PropsSI('D', 'T', T_in, 'P', p, fluid)
    m = rho_in*U_inf*N_T*S_T
    ARF = 0.1
    T_out = T_out0
    
    for i in range(N_max):
        # Save the old T_out from last iteration round
        T_out_old = T_out 
        # Calculate new values for the variables
        h = h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid) 
        Q = N*h*sp.pi*D*dT_log(T_in, T_out, T_s)
        # Calculate new T_out'''
        c_p = PropsSI('C', 'T', (T_in+T_out)/2, 'P', p, fluid)
        T_out_new = T_in + Q/(m*c_p)
        # Calculate error and stop the for-loop if within tolerance
        error = abs( T_out_new - T_out_old )
        if error < 0.01: # Limit for difference between two iterations (K)
            print('Iterations required:', i)
            T_out = T_out_new
            break
        # Use underrelaxation for next iteration round
        T_out = ARF*T_out_new + (1-ARF)*T_out_old
        
    return h, Q, T_out

#---------------------  MAIN CODE --------------------------
# Initial values for calculations
D = 0.635
S_T = 2.5*D
N_T = 24
N_L = sp.linspace(20, 300, 50)
N = N_T*N_L
T_in = 900
T_s = 600
U_inf = 5
p = 1e5
fluid = 'Air'

N_max = 500 # Maximum number of iterations
T_out0 = (T_s+T_in)/2 # Initial guess for outlet temperature
h, Q, T_out = sp.zeros(len(N_L)), sp.zeros(len(N_L)), sp.zeros(len(N_L))

for i in range(len(N_L)):
    h[i], Q[i], T_out[i] = iterate_solution(T_out0, N_max, N[i], S_T, U_inf, 
                                              D, T_s, T_in, p, fluid)

plt.figure(0)      
plt.plot(N_L, T_out)
plt.xlabel("Number of tube columns #")
plt.ylabel("Fluid outlet temperature (K)")
plt.grid('On')

plt.figure(1)
plt.plot(N_L, abs(Q)/1e6)
plt.xlabel("Number of tube columns #")
plt.ylabel("Heat flux per meter of pipe (MW/m)")
plt.grid('On')


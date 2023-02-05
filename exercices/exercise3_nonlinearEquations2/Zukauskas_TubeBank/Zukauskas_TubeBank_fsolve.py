# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 11:50:58 2019

@author: Niko Niemelä
"""
import scipy as sp
from scipy.optimize import fsolve, brentq
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI

#----  Zukauskas correlation for inline tube banks   ----

def h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid):

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

#--------------------- dT_log --------------------------

def dT_log(T_in, T_out, T_s):
    
    dT_log = ( (T_s-T_in)-(T_s-T_out) ) / sp.log( (T_s-T_in)/(T_s-T_out) )
    return dT_log

#---- Heat transfer rate per meter of pipe length, Q [W/m] ----

def Q_bundle(N, S_T, U_inf, D, T_in, T_out, T_s, p, fluid):

    h = h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid) 
    Q_bundle = N*h*sp.pi*D*dT_log(T_in, T_out, T_s)
    return Q_bundle

#---------------------  Iterate solution   --------------------------
    
def calculate_solution(T_out0, N, S_T, U_inf, D, T_in, T_s, p, fluid):
    
    rho_in = PropsSI('D', 'T', T_in, 'P', p, fluid)
    m = rho_in*U_inf*N_T*S_T

    def to_zero(T_out):
        T_m = (T_in+T_out)/2
        c_p = PropsSI('C', 'T', T_m, 'P', p, fluid)
        Q = Q_bundle(N, S_T, U_inf, D, T_in, T_out, T_s, p, fluid)
        return Q - m*c_p*(T_out-T_in)

    T_out = brentq(to_zero, T_s+0.1, T_in-0.1)
    Q = Q_bundle(N, S_T, U_inf, D, T_in, T_out, T_s, p, fluid)
    h = h_m(S_T, U_inf, D, T_in, T_out, T_s, p, fluid)
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

# fsolve needs a good initial guess for each N_L. Good guess can be obtained
# using the exponential formula with decent estimate for h.
rho_in = PropsSI('D', 'T', T_in, 'P', p, fluid)
m = rho_in*U_inf*N_T*S_T
h_0 = h_m(S_T, U_inf, D, T_in, (T_in+T_s)/2, T_s, p, fluid)
c_p = PropsSI('C', 'T', (T_in+T_s)/2 , 'P', p, fluid)
# Initial guess for outlet temperature for each N_L
T_out0 = (T_in - T_s)*sp.exp( (-N*h_0*sp.pi*D) / (m*c_p) ) + T_s

h, Q, T_out = sp.zeros(len(N_L)), sp.zeros(len(N_L)), sp.zeros(len(N_L))

for i in range(len(N_L)):
    h[i], Q[i], T_out[i] = calculate_solution(T_out0[i], N[i], S_T, U_inf, D, T_in, T_s, p, fluid)

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


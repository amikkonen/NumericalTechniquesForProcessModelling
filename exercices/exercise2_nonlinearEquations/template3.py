# -*- coding: utf-8 -*-
import numpy as np
from scipy import optimize
from CoolProp.CoolProp import PropsSI


C2K = 273.15
# emissivity 
eps = 0.8
# Stefan-Boltman constant
sigma = 5.67036713e-8
# Graaviational acceleration
g = 9.81
# P.S Python3 acepts any unicode character, you could write 
#σ = 5.67036713e-8

# Parameters
L_vertical = 10e-2
L_horizontal = 5e-2

Q1 = 20 # Watts
Q2 = 40 # Watts
T_inf = 20+C2K
F = 0.3
fluid = "air"
p = 1e5

A = L_vertical*L_horizontal

def hVertical(T):
    """
    Natural convection for vertical plate.
    
    Correlation grapped from:
        https://en.wikipedia.org/wiki/Natural_convection#Natural_convection_from_a_vertical_plate
    Applicapility or primary sources have not been checked! 
    Do not use in real applications without chekking! 
    
    For convinience, we use CoolProp here for fluid properties. In real 
    applications it would be adviceble to find same faster solution. CoolProps 
    is slow and deep inside iterative loops performace matters.
    """
    # density
    rho = PropsSI("D", "T", T, "P", p, fluid)
    # dynamic viscosity, directly available in CoolProps
    mu  = PropsSI("V", "T", T, "P", p, fluid)
    # Kinematic viscosity
    nu  = mu/rho
    # Heat conductivity
    kf  = PropsSI("conductivity", "T", T, "P", p, fluid)
    # Coefficient of thermal expansion (equal to approximately 1/T, for ideal gases)
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", p, fluid)
    
    # Grashof number, a typical dimensionless number for natural convection
    Gr = (beta*g*L_vertical**3*(T-T_inf)) / (nu**2)
    # Nusslet number
    Nu = 0.478*Gr**0.25
    # Heat transfer coefficient
    h = Nu * kf / L_vertical
    
    return h
    
def energy_balance(T):
    """Returns energy inbalance for object 1 and 2. Make this zero.
    
    """
    
    # Calculate convective heat transfer coefficients
    # CoolProps does not accept arrays, therefore we need 2 separate calls
    h1 = hVertical(T[0])
    h2 = hVertical(T[1])
    
    # Calculate energy balance
    f = np.zeros(2)
#    f[0] = (
#                          # Object radiates out
#                          # The other object radiates in
#                          # Enviroment radiates in
#                          # Convective heat transfer out
#                          # Electric heating in
#           )
#    f[1] = (
#           )
    return f
    
# Solve T with fsolve
T_initial = np.array([50,50])+C2K
#T = write_code_here

# Output
#print(T-C2K)


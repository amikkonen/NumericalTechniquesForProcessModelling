# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
from CoolProp.CoolProp import PropsSI

# Celsius to Kelvin
C2K = 273.15
# Stefan-Boltman constant
sigma = 5.67036713e-8
# Graaviational acceleration
g = 9.81
# P.S Python3 acepts any unicode character, you could write 
#σ = 5.67036713e-8



def haaland(Re, eps_over_d):
    """ Friction factor, Haaland equation
    """
    f  = (-1.8*np.log10((eps_over_d/3.7)**1.11+6.9/Re))**-2
    return f

def heat_transfer_coefficient_natural_convection_vertical_plate(T, T_inf, 
                                                                L, p, fluid):
    """
    Natural convection for vertical plate.
    
    Correlation grapped from:
        https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Natural%20Convection.pdf
    Applicapility or primary sources have not been checked! 
    Do not use in real applications without chekking! 
    
    For convinience, we use CoolProp here for fluid properties. In real 
    applications it would be adviceble to find a faster solution. CoolProps 
    is slow and deep inside iterative loops performace matters.
    """
    T_m = (T+T_inf)/2
    
    # density
    rho = PropsSI("D", "T", T_m, "P", p, fluid)
    # dynamic viscosity, directly available in CoolProps
    mu  = PropsSI("V", "T", T_m, "P", p, fluid)
    # Kinematic viscosity
    nu  = mu/rho
    # Heat conductivity
    kf  = PropsSI("conductivity", "T", T_m, "P", p, fluid)
    # Coefficient of thermal expansion (equal to approximately 1/T, for ideal gases)
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_m, "P", p, fluid)
    # Prandtl number
    Pr  = PropsSI("Prandtl", "T", T_m, "P", p, fluid)
    
    # Grashof number, a typical dimensionless number for natural convection
    Gr = (beta*g*L**3*(T-T_inf)) / (nu**2)
    
    # Rayleigh
    Ra = Gr*Pr
    # Nusselt number
    if Ra < 1e4 or Ra > 1e13:
        raise ValueError("Ra < 1e4 or Ra > 1e13", "Ra=", Ra)
    elif Ra < 1e9:
        Nu = 0.59*Ra**0.25
    else:
        Nu = 0.1*Ra**(1/3)
    
    # Heat transfer coefficient
    h = Nu * kf / L
    
    return h
    
def heat_transfer_coefficient_forced_convection_duct(V, D, T, T_in, p, 
                                                     eps_over_d, fluid):
    """Gnielinski's correlation for turbulent flow in tubes
    
    """
    T_m = (T + T_in)/2
    
    
    # Fluid properties    
    # density
    rho = PropsSI("D", "T", T_m, "P", p, fluid)
    # dynamic viscosity, directly available in CoolProps
    mu  = PropsSI("V", "T", T_m, "P", p, fluid)
    # Kinematic viscosity
    nu  = mu/rho
    # Prandtl number
    Pr  = PropsSI("Prandtl", "T", T_m, "P", p, fluid)
    # Heat conductivity
    kf  = PropsSI("conductivity", "T", T_m, "P", p, fluid)
    
    # Reynolds number
    Re = V*D/nu
    
    # Friction factor, Haaland equation
    f  = haaland(Re, eps_over_d)
    # Nusselt number
    Nu = ((f/8)*(Re-1000)*Pr) / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
    # Heat transfer coefficient
    h  = Nu*kf/D

    return h

def heat_transfer_coefficient_vertical_infinite_cavity(T, p, L, fluid):
    """Heat transfer coefficient in infinite vertical cavity.
    
    T - arrays len 2, the plate temperatures
    
    Building Physics - Heat, Air and Moisture, Hugo Hens, Eq.1.59, pp.57
    
    """
    
    T_m = T.mean()
    
    # density
    rho = PropsSI("D", "T", T_m, "P", p, fluid)
    # dynamic viscosity, directly available in CoolProps
    mu  = PropsSI("V", "T", T_m, "P", p, fluid)
    # Kinematic viscosity
    nu  = mu/rho
    # Coefficient of thermal expansion (equal to approximately 1/T, for ideal gases)
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_m, "P", p, fluid)
    # Prandtl number
    Pr  = PropsSI("Prandtl", "T", T_m, "P", p, fluid)
    # Heat conductivity
    kf  = PropsSI("conductivity", "T", T_m, "P", p, fluid)
    
    Ra = (beta*g*L**3*(np.absolute(T[0]-T[1]) )) / (nu**2) * Pr
    Nu = max(1, 1+(0.024*Ra**1.39) / (Ra+10100) )


    # Heat transfer coefficient
    h  = Nu*kf/L

    return h


#############################################################################
# MAIN
#############################################################################

# Parameters

# emissivity 
eps = 0.9
# Surface roughness
eps_over_d = 0
# Velocity inside duct
V = 2
# Hydraulic diameter 
D = 0.3
# Duct hight
L_v = 2
# Cavity size
L_c = 1e-2
# Temperature inside duct 
T_in  = 800+C2K
# Temperature outside
T_inf = 20+C2K
# Reference pressure
p = 1e5
# We use air as the fluid
fluid = "air"

def energy_balance(T):
    """Returns energy inbalance for plates 1 and 2. Zero when correct.
    
    """
    # Heat transfer coefficient inside duct
    h_i0 = heat_transfer_coefficient_forced_convection_duct(
            V, D, T[0], T_in, p, eps_over_d, fluid)
    # Heat transfer coefficient inside cavity
    h_01 = heat_transfer_coefficient_vertical_infinite_cavity(
            T, p, L_c, fluid)
    # Heat transfer coefficient outside
    h_1inf = heat_transfer_coefficient_natural_convection_vertical_plate(
            T[1], T_inf, L_v, p, fluid)
    
    # Calculate energy balance
    f = np.zeros(2)
    f[0] = (
                h_i0*(T_in - T[0] )                             # convection in duct
               - h_01 *(T[0] - T[1] )                # convection in cavity
               - eps*sigma*(T[0]**4 - T[1]**4)                             # radiation in cavity
          )
    f[1] = (
               h_01*(T[0]-T[1])             # convection in cavity
             - h_1inf*(T[1]-T_inf)          # convection outside
             + eps*sigma*(T[0]**4-T[1]**4)  # radiation in cavity
             - eps*sigma*(T[1]**4-T_inf**4) # radiation outside
           )
    # Units W/m^2
    return f



########################################
# (a) Try to guess the temperatures
########################################

# T_initial = np.array([260, 160])+C2K
# print("Energy inbalance", energy_balance(T_initial))


########################################
# (b) Use root
########################################

T_initial = np.array([300, 200])+C2K
res = optimize.root(energy_balance, T_initial)
#    print(res)
assert(res.success)
T = res.x
print("Temperatures (C)", (T-C2K).round(1))

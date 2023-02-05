# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from CoolProp.CoolProp import PropsSI

# Parameters
# Remember to use SI units
# You can use a short hand for small and large numbers
# 1e2 == 100, 1e5 == 100000, 4e-3 == 0.004 etc...
C2K = 273.15
D = 10e-3
Q = 15/60/1000 # l/min -> SI
T_in = 30 + C2K
T_w = 120 + C2K
T_out = 90 + C2K
p_out = 2e5

# Mean temperature
T_m = (T_in + T_out)/2

# Fluid properties    
# http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
fluid = "water"
#fluid = "air"
# density
rho = PropsSI("D", "T", T_m, "P", p_out, fluid)
# dynamic viscosity, directly available in CoolProps
mu  = PropsSI("V", "T", T_m, "P", p_out, fluid)
# Kinematic viscosity
nu  = mu/rho
# Prandtl number
Pr  = PropsSI("Prandtl", "T", T_m, "P", p_out, fluid)
# Heat conductivity
kf  = PropsSI("conductivity", "T", T_m, "P", p_out, fluid)
# Heat capasity
cp  = PropsSI("C", "T", T_m, "P", p_out, fluid)

# Mass flow rate, pipe cross section area, and npeed
m  = Q*rho
A  = 0.25*np.pi*D**2
V  = Q/A

# Reynolds number
Re = V*D/nu
# Friction factor, Haaland equation, eq. 3
f  = (-1.8*np.log10(6.9/Re))**-2
# Nusselt number, eq. 2
Nu = ((f/8)*(Re-1000)*Pr) / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
# Heat transfer coefficient
h  = Nu*kf/D
# Pipe length, eq. 1
L = -np.log((T_w -T_out) / (T_w-T_in)) / ( h*2*np.pi*(D/2) / (m*cp) )
# Pressure difference and pressure at inlet, eq 4
dp = 0.5*rho*V**2*L/D*f
p_in = p_out + dp

print("Re", Re)
print("L", L)
print("p_in (bar)", p_in*1e-5) # For output, convert Pa -> bar

# Pipe length for plotting
# https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.linnpace.html
L_ar = np.linspace(1e-5, L)
# Calculate temperatures at different locations
T_fluid = T_w - (T_w-T_in) * np.exp( -h*2*np.pi*(D/2)*L_ar / (m*cp) )

# Plotting
plt.plot(L_ar, T_fluid-C2K)  # For output, convert K -> C
plt.xlabel("x (m)")
plt.ylabel("T (C)")
plt.grid(True)

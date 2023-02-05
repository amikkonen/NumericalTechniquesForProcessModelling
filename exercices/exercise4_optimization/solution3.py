# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
from CoolProp.CoolProp import PropsSI

# Celsius to Kelvin
C2K = 273.15

p_out = 1e5
T_in  = 30+C2K
T_out = 80+C2K
T_w   = 100 +C2K
fluid = "water"

thickness = 2e-3
rho_pipe = 8960 # Copper denity
m = 0.5
dp_max = 1e3

# Mean temperature
T_m = (T_in + T_out)/2

# Fluid properties    
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


# the argument final=False defaults to False and can be ignored

def pipe(D, final=False):
    
    # Volumetric flow rate, area, and speed
    A  = 0.25*np.pi*D**2
    Q  = m/rho 
    V  = Q/A
    
    # Reynolds number
    Re = V*D/nu
    # Friction factor, Haaland equation
    f  = (-1.8*np.log10(6.9/Re))**-2
    # Nusselt number
    Nu = ((f/8)*(Re-1000)*Pr) / (1 + 12.7*(f/8)**0.5*(Pr**(2/3)-1))
    # Heat transfer coefficient
    h  = Nu*kf/D
    # Pipe length
    L = -np.log((T_w -T_out) / (T_w-T_in)) / ( h*2*np.pi*(D/2) / (m*cp) )
    # Pressure difference and pressure at inlet
    dp = 0.5*rho*V**2*L/D*f
    # Pipe mass
    m_pipe = 0.25*np.pi*((D+2*thickness)**2-D**2)*L*rho_pipe

    # If final, print Reynolds number
    if final:
        print("Re", Re)

    # Return pipe mass, pressure loss and pipe length
    return m_pipe, dp, L

def pipe_mass_with_penalty(D):
    """The function to be optimized. Define a suitable (large) penalty to 
    change the location of minumum to fit the constrain.
    
    Returns pipe mass as a function of dipe diameter.
    """
    m_pipe, dp, L = pipe(D)
    
    if dp > dp_max:
        m_pipe += 100 
    return m_pipe



res = optimize.minimize_scalar(pipe_mass_with_penalty, bounds=[1e-3, 10e-2], 
                               method="bounded")
# Print optimal solution and
# check for reasonable Reynolds number
print("mass (kg), dp (Pa), L (m)", pipe(res.x, final=True))
print()

# Print results
print(res)

###########################################################################
# Plotting
###########################################################################
def pipe_mass(D):
    m_pipe, dp, L = pipe(D)
    return m_pipe

def pressure_loss(D):
    m_pipe, dp, L = pipe(D)
    return dp

d = np.linspace(3e-2, 5e-2)

# pipe mass
plt.plot(d*1e2, pipe_mass(d), label="pipe mass")
plt.plot(res.x*1e2, pipe_mass(res.x), "dr", label="optimal")

plt.xlabel("d (cm)")
plt.ylabel("m (kg)")
plt.legend(frameon=False)

# Pressure loss
twinx = plt.twinx()
twinx.plot(d*1e2, pressure_loss(d), "--", label="pressure loss")

# "Graphical solution"
twinx.axhline(pressure_loss(res.x), color="g", label="max pressure loss")
twinx.axvline(res.x*1e2, color="y",label="\"graphical solution\"")
twinx.set_ylabel("dp (Pa)")
twinx.legend(frameon=False)
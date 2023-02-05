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


## Fluid properties    
#
## the argument final=False defaults to False and can be ignored
#def pipe(D, final=False):
#
#
#    # If final, print Reynolds number
#    if final:
#        print("Re", Re)
#
#    # Return pipe mass, pressure loss and pipe length
#    return m_pipe, dp, L

#def pipe_mass_with_penalty(D):
#    """The function to be optimized. Define a suitable (large) penalty to 
#    change the location of minumum to fit the constrain.
#    
#    Returns pipe mass as a function of dipe diameter.
#    """
#    m_pipe, dp, L = pipe(D)
#    
#    if 
#        m_pipe +=  
#    return m_pipe

# Testing
#print(pipe(5e-2))

# Optimization
#res = 

# Print optimal solution and
# check for reasonable Reynolds number
#print("mass (kg), dp (Pa), L (m)", pipe(res.x, final=True))
#print()
#
## Print results
#print(res)


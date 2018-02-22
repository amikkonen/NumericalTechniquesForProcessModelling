#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

KEB-45250
Numerical Techniques for Process Modelling,
Prosessien numeerinen mallinnus

Exercice 6

Emphesis on simplisity. Low performance in larger systems.

Created on 2018 Feb 9 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""
    
##############################################################################

# Lenght (m)
L = 0.5
# Pipe diameter (m)
d = 0.02
# Inlet velocity (m/s)
V = 5
# Viscosity
nu = 0.01
# Density (kg/m3)
rho = 1

# Reynolds number
Re = V*d/nu
print("Re", Re)

# Pressure drop
f = 64/Re
dp = 0.5*f*rho*V**2*L/d 

# Print result
print("dp", dp)


